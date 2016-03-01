function s = PriorAnalysis(s)

% PriorAnalysis
%
% Analyzes the priors
%
% See also:
%
% ...........................................................................
% 
% Created: January 28, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

Action = 'PriorAnalysis';

% % check if model already prepared
% if isfield(s.Status,Action) && s.Status.(Action), return, end

fprintf('\n*** Analyzing DSGE Prior distribution\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% Options

s = CheckOptions(s);

if isfield(s.Options,'Prior')
    op = s.Options.Prior;
else
    op = struct; 
end
if ~isfield(op,'nDraws'),op.nDraws = 1000; end
if ~isfield(op,'Percentiles')
    op.Percentiles = [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99];
end

s.FileName.PriorDraw = [s.Spec,'_PriorDraw'];
% s.FileName.PriorSample = [s.Spec,'_PriorSample'];

s.Report.Prior = sprintf('%s_Report_Prior',s.Spec);
ReportTitle = sprintf('Prior Analysis:\\\\%s',s.Spec);

s.Options.Prior = op;

opTable = s.Options.ParTable;

%% -------------------------------------------------------------------

%% Prepare variables
nprc = length(op.Percentiles);
PrcList = cell(nprc,1);
for jprc=1:nprc
    PrcList{jprc} = sprintf('PriorPrc%03.0f',1000*op.Percentiles(jprc));
end


%% Analyze Parameters
p = s.Param;
np = s.n.Param;
p.PriorMode = nan(np,1);
p.PriorParams = nan(np,2);
for jprc=1:nprc
    p.(PrcList{jprc}) = nan(np,1);
end
p.PriorLPdfCmd = cell(np,1);
p.PriorPdfCmd = cell(np,1);
p.PriorRndCmd = cell(np,1);
for j=1:np
    if strcmp(p.PriorDist{j},'N')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        p.PriorMode(j) = pmean;
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = norminv(op.Percentiles(jprc),pmean,psd);
        end
        p.PriorParams(j,:) = [pmean,psd];
        p.PriorLPdfCmd{j} = sprintf('log(normpdf(%s,%.16f,%.16f))',...
                                    p.Names{j},pmean,psd);
        p.PriorPdfCmd{j} = sprintf('normpdf(%%.16f,%.16f,%.16f)',pmean,psd);
        p.PriorRndCmd{j} = sprintf('normrnd(%.16f,%.16f)',pmean,psd);
    elseif strcmp(p.PriorDist{j},'TN')
        % Assume x>=0
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        a = -pmean/psd;
        acdf = normcdf(a,0,1);
        aZ = 1-normcdf(a,0,1);
        alambda = normpdf(a,0,1)/aZ;
        adelta = alambda*(alambda-a);
        p.PriorMean(j) = pmean + psd*alambda;
        p.PriorSD(j) = psd*(1-adelta)^(1/2);
        p.PriorMode(j) = max(0,pmean);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = norminv(...
                op.Percentiles(jprc)*aZ+acdf,pmean,psd);
        end
        p.PriorParams(j,:) = [pmean,psd];
        p.PriorLPdfCmd{j} = sprintf(...
            ['log((%1$s>=0)'...
             '*normpdf((%1$s-%2$.16f)/%3$.16f,0,1)/%3$.16f/%4$.16f)'],...
            p.Names{j},pmean,psd,aZ);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>=0)',...
             '*normpdf((%%1$.16f-%1$.16f)/%2$.16f,0,1)/%2$.16f/%3$.16f'],...
            pmean,psd,aZ);
        p.PriorRndCmd{j} = sprintf('norminv(rand*%.16f+%.16f,%.16f,%.16f)',...
                                   aZ,acdf,pmean,psd);
    elseif strcmp(p.PriorDist{j},'B')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        a = pmean*(pmean-pmean^2-psd^2)/psd^2;
        b = a*(1/pmean-1);
        p.PriorMode(j) = min(max(0,(a-1)/(a+b-2)),1);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = betainv(op.Percentiles(jprc),a,b);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf('log(betapdf(%s,%.16f,%.16f))',...
                                    p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf('betapdf(%%.16f,%.16f,%.16f)',a,b);
        p.PriorRndCmd{j} = sprintf('betarnd(%.16f,%.16f)',a,b);
    elseif strcmp(p.PriorDist{j},'G')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        a = (pmean/psd)^2;
        b = pmean/a;
        if a>=1
            p.PriorMode(j) = (a-1)*b;
        else
            p.PriorMode(j) = NaN;
        end
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(op.Percentiles(jprc),a,b);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf('log(gampdf(%s,%.16f,%.16f))',...
                                    p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf('gampdf(%%.16f,%.16f,%.16f)',a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)',a,b);
    elseif strcmp(p.PriorDist{j},'IG1')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        if psd==inf
            a = 1;
        else
            fname = sprintf('igamsolve%.0f',cputime*1e10);
            fid=fopen([fname,'.m'],'wt');
            fprintf(fid,'function f=%s(x)\n',fname);
            fprintf(fid,'pmean = %.16f;\n',pmean);
            fprintf(fid,'pvar = %.16f;\n',psd^2);
            fprintf(fid,'for j=1:length(x)\n');
            fprintf(fid,'    a = x(j);\n');
            fprintf(fid,...
                    ['    f(j) = 1/(a-1)*(pmean*gamma(a)/gamma(a-1/2))^2',...
                     '-pmean^2-pvar;\n']);
            fprintf(fid,'end\n');
            fclose(fid);
            [a,rc] = csolvevb(fname,5,[],1e-10,1000);
            if rc~=0, 
                error('Search for iGam parameters failed, rc=%.0f',rc), 
            end
            delete([fname,'.m'])
        end
        b = (gamma(a-1/2)/pmean/gamma(a))^2;
        p.PriorMode(j) = (1/b/(a+1/2))^(1/2);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(1-op.Percentiles(jprc),a,b)^(-1/2);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-2),%2$.16f,%3$.16f)*2/%1$s^3))',...
            p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-2),%1$.16f,%2$.16f)*2/%%1$.16f^3)'],a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)^(-1/2)',a,b);
    elseif strcmp(p.PriorDist{j},'IG2')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        if psd==inf
            a = 2;
        else
            a = 2+pmean^2/psd^2;
        end
        b = 1/pmean/(a-1);
        p.PriorMode(j) = 1/b/(a+1);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(1-op.Percentiles(jprc),a,b)^(-1);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',...
            p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)'],a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)^(-1)',a,b);
    elseif strcmp(p.PriorDist{j},'C')
        pmean = p.PriorMean(j);
        psd = 0;
        p.PriorSD(j) = 0;
        p.PriorMode(j) = pmean;
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = pmean;
        end
        p.PriorParams(j,:) = [pmean,psd];
        p.PriorLPdfCmd{j} = sprintf('log(1*(%s==%.16f))',p.Names{j},pmean);
        p.PriorPdfCmd{j} = sprintf('1*(%%.16f==%.16f)',pmean);
        p.PriorRndCmd{j} = sprintf('%.16f',pmean);
    end
end
s.Param = p;
Prior.UnconstrainedParam = p;

%% display results on screen
fprintf('\nPrior (Unconstrained):')
fprintf('\n----------------------\n')
namelength = [cellfun('length',p.Names)];
namelengthmax = max(namelength);
DispList = {'','Names';
            'Dist','PriorDist';
            '  Mode','PriorMode';
            '  Mean','PriorMean';
            '   SD','PriorSD';
            '   5%','PriorPrc050';
            ' Median','PriorPrc500';
            '   95%','PriorPrc950';
           }';
nc = size(DispList,2);
fprintf(['%-',int2str(namelengthmax),'s'],DispList{1,1});
fprintf('  %-4s',DispList{1,2});
for jc=3:nc
    fprintf('  %-8s',DispList{1,jc});
end
fprintf('\n');
for j=1:np
    fprintf(['%',int2str(namelengthmax),'s'],p.(DispList{2,1}){j});
    fprintf('  %4s',p.(DispList{2,2}){j});
    for jc=3:nc
        fprintf('  %8.4f',p.(DispList{2,jc})(j));
    end
    fprintf('\n');
end
fprintf('\n');


%% Generate prior draws

fprintf('Generating function to draw from prior...\n')

fid = fopen([s.FileName.PriorDraw,'.m'],'wt');
fprintf(fid,'function x=%s(nDraws)\n\n',s.FileName.PriorDraw);
fprintf(fid,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs\n\n',clock);
fprintf(fid,'if ~exist(''nDraws'',''var''), nDraws = 1; end\n');
fprintf(fid,'x = zeros(%.0f,nDraws);\n',np);
fprintf(fid,'for jd=1:nDraws\n');
for j=1:np
    fprintf(fid,'    x(%.0f,jd) = %s;\n',j,p.PriorRndCmd{j});
end
fprintf(fid,'end\n');
fclose(fid);

fprintf('\nPrior Sample:')
fprintf('\n-------------\n\n')
Prior.nDraws = op.nDraws;
xd = zeros(np,Prior.nDraws);
xdAux = zeros(s.n.AuxParam,Prior.nDraws);
xj = zeros(np,1);
BadDraws = false(1,Prior.nDraws);
xd = feval(s.FileName.PriorDraw,Prior.nDraws);
fn = s.FileName.Mats;
nAux = s.n.AuxParam;
ListAux = s.AuxParam.Names;
parfor jd=1:Prior.nDraws
    Matsj = feval(fn,xd(:,jd));
    BadDraws(jd) = ~all(Matsj.REE.eu==1) || Matsj.KF.sig00rc~=0;
    for jp=1:nAux
        xdAux(jp,jd) = Matsj.AuxParam.(ListAux{jp});
    end
end
clear fn
Prior.nBadDraws = sum(BadDraws);
Prior.FractionBadDraws = Prior.nBadDraws/Prior.nDraws;
Prior.LogTruncationCorrection = -log(1-Prior.FractionBadDraws);
xd(:,BadDraws) = [];
Prior.nDraws = size(xd,2);
fprintf('Number of accepted draws: %.0f\n',Prior.nDraws);
fprintf(...
    'Percent of rejected draws: %.2f%%\n',...
    Prior.FractionBadDraws*100);
fprintf('log-prior correction: %.6f\n',Prior.LogTruncationCorrection);

p.PriorMean = mean(xd,2);
p.PriorSD = std(xd,0,2);
for jprc=1:nprc
    p.(PrcList{jprc}) = prctile(xd,100*op.Percentiles(jprc),2);
end
Prior.Param = p;
s.Param = p;

xdAux(:,BadDraws) = [];
pAux = s.AuxParam;
pAux.PriorMean = mean(xdAux,2);
pAux.PriorSD = std(xdAux,0,2);
for jprc=1:nprc
    pAux.(PrcList{jprc}) = prctile(xdAux,100*op.Percentiles(jprc),2);
end
s.AuxParam = pAux;
Prior.AuxParam = pAux;

s.Prior = Prior;

pList = {'Param','AuxParam'};
DispList = {'  Mean','PriorMean';
            '   SD','PriorSD';
            '   5%','PriorPrc050';
            ' Median','PriorPrc500';
            '   95%','PriorPrc950';
           }';
nc = size(DispList,2);
for jP=1:length(pList)
    Pj = pList{jP};
    psj = Prior.(Pj);
    namelength = [cellfun('length',s.(Pj).Names)];
    namelengthmax = max(namelength);
    fprintf(['\n%-',int2str(namelengthmax),'s'],'');
    for jc=1:nc
        fprintf('  %-8s',DispList{1,jc});
    end
    fprintf('\n');
    for jp=1:s.n.(Pj)
        fprintf(['%',int2str(namelengthmax),'s'],s.(Pj).Names{jp});
        for jc=1:nc
            fprintf('  %8.4f',Prior.(Pj).(DispList{2,jc})(jp));
        end
        fprintf('\n');
    end
end
fprintf('\n');
% save(s.FileName.PriorSample,'xd','xdAux')


%% Make Prior Report

fprintf('Making report: %s\n',s.Report.Prior);
fid = vcCreateTex(s.Report.Prior,ReportTitle);
fprintf(fid,'\\newpage \n');
fprintf(fid,'\\section{Parameters}\n');
str = [' & %.',int2str(opTable.Precision),'f'];
TableBreaks = opTable.MaxRows:opTable.MaxRows:s.n.Param;
if ~ismember(s.n.Param,TableBreaks), TableBreaks(end+1) = s.n.Param; end
idxPar = 0;
nBreaks = length(TableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    if opTable.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{lcccccccccccc} \n');
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{7}{c}{Unconstrained Prior} ');
    fprintf(fid,'& & \\multicolumn{4}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mean & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'$%s$',s.Param.PrettyNames{jr});
        fprintf(fid,' & %s', s.Param.PriorDist{jr});
        fprintf(fid,str,Prior.UnconstrainedParam.PriorMode(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorMean(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorSD(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc050(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc500(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc950(jr));
        fprintf(fid,' &');
        fprintf(fid,str,Prior.Param.PriorMean(jr));
        fprintf(fid,str,Prior.Param.PriorPrc050(jr));
        fprintf(fid,str,Prior.Param.PriorPrc500(jr));
        fprintf(fid,str,Prior.Param.PriorPrc950(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,opTable.Lines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\section{Auxiliary Parameters}\n');
str = [' & %.',int2str(opTable.Precision),'f'];
TableBreaks = opTable.MaxRows:opTable.MaxRows:s.n.AuxParam;
if ~ismember(s.n.Param,TableBreaks), TableBreaks(end+1) = s.n.AuxParam; end
idxPar = 0;
nBreaks = length(TableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Auxiliary Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{lccccccccccccc} \n');
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Mean & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'$%s$',s.AuxParam.PrettyNames{jr});
        fprintf(fid,str,Prior.AuxParam.PriorMean(jr));
        fprintf(fid,str,Prior.AuxParam.PriorPrc050(jr));
        fprintf(fid,str,Prior.AuxParam.PriorPrc500(jr));
        fprintf(fid,str,Prior.AuxParam.PriorPrc950(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,opTable.Lines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(s.Report.Prior)


%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);
fprintf('\n%s %s\n\n',Action,vctoc([],s.TimeElapsed.(Action)))

%% -------------------------------------------------------------------

