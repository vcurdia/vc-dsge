function obj = AnalyzePrior(obj)

% AnalyzePrior
%
% Analyzes the priors
%
% See also:
% DSGE, SetupMyDSGE
%
% ...........................................................................
% 
% Created: November 10, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

action = 'AnalyzePrior';
obj = obj.TrackTime(action,1);

fprintf('\n*** Analyzing DSGE Prior distribution\n')

%% Options

if ~isempty(obj.Prior)
    Prior = obj.Prior;
else
    Prior = struct; 
end

obj.FileName.DrawPrior = [obj.Name,'_DrawPrior'];

ReportFileName = sprintf('%s_Report_Prior',obj.Name);
ReportTitle = sprintf('Prior Analysis:\\\\%s',obj.Name);

%% -------------------------------------------------------------------

%% Prepare variables
nPrc = length(obj.ParamPercentiles);
PrcList = cell(nPrc,1);
for jPrc=1:nPrc
    PrcList{jPrc} = sprintf('PriorPrc%03.0f',1000*obj.ParamPercentiles(jPrc));
end

%% Analyze Parameters
p = obj.Param;
np = p.N;
p.PriorMode = nan(np,1);
p.PriorParams = nan(np,2);
for jPrc=1:nPrc
    p.(PrcList{jPrc}) = nan(np,1);
end
p.PriorLPdfCmd = cell(np,1);
p.PriorPdfCmd = cell(np,1);
p.PriorRndCmd = cell(np,1);
pOptions = optimoptions(@fsolve);
pOptions.Display = 'off';
for j=1:np
    if strcmp(p.PriorDist{j},'N')
        pmean = p.PriorMean(j);
        psd = p.PriorSD(j);
        p.PriorMode(j) = pmean;
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = norminv(obj.ParamPercentiles(jPrc),pmean,psd);
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
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = norminv(...
                op.Percentiles(jPrc)*aZ+acdf,pmean,psd);
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
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = betainv(obj.ParamPercentiles(jPrc),a,b);
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
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = gaminv(obj.ParamPercentiles(jPrc),a,b);
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
%             [a,rc] = csolvevb(@(x)igamsolve(x,pmean,psd),5,[],1e-10,1000);
            [a,~,rc] = fsolve(@(x)igamsolve(x,pmean,psd),5,pOptions);
%             if rc~=0, 
            if rc~=1, 
                error('Search for iGam parameters failed, rc=%.0f',rc), 
            end
        end
        b = (gamma(a-1/2)/pmean/gamma(a))^2;
        p.PriorMode(j) = (1/b/(a+1/2))^(1/2);
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = gaminv(1-obj.ParamPercentiles(jPrc),a,b)^(-1/2);
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
        for jPrc=1:nPrc
            p.(PrcList{jPrc})(j) = gaminv(1-obj.ParamPercentiles(jPrc),a,b)^(-1);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',...
            p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)'],a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)^(-1)',a,b);
    end
end
obj.Param = p;
Prior.UnconstrainedParam = p;

%% display results on screen
fprintf('\nPrior (Unconstrained):')
fprintf('\n----------------------\n')
namelength = [cellfun('length',p.Names)];
namelengthmax = max(namelength);
DispList = {'','Names';
            'Dist','PriorDist';
            '    Mode','PriorMode';
            '    Mean','PriorMean';
            '      SD','PriorSD';
            '      5%','PriorPrc050';
            '  Median','PriorPrc500';
            '     95%','PriorPrc950';
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

fid = fopen([obj.FileName.DrawPrior,'.m'],'wt');
fprintf(fid,'function x=%s(nDraws)\n\n',obj.FileName.DrawPrior);
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
xd = zeros(np,obj.PriorNDraws);
nAux = obj.AuxParam.N;
xdAux = zeros(nAux,obj.PriorNDraws);
xj = zeros(np,1);
BadDraws = false(1,obj.PriorNDraws);
xd = obj.DrawPrior(obj.PriorNDraws);
fh = @(x)obj.Mats(x);
AuxNames = obj.AuxParam.Names;
% Matsd = cell(obj.PriorNDraws);
parfor jd=1:obj.PriorNDraws
%     jd
    Matsj = fh(xd(:,jd));
    BadDraws(jd) = ~Matsj.Status==1;
    for jp=1:nAux
        xdAux(jp,jd) = Matsj.AuxParam.(AuxNames{jp});
    end
%     Matsd{jd} = Matsj;
end
% keyboard
Prior.NDraws = obj.PriorNDraws;
Prior.NBadDraws = sum(BadDraws);
Prior.FractionBadDraws = Prior.NBadDraws/obj.PriorNDraws;
Prior.LogTruncationCorrection = -log(1-Prior.FractionBadDraws);
xd(:,BadDraws) = [];
Prior.NDraws = size(xd,2);
fprintf('Number of accepted draws: %.0f\n',Prior.NDraws);
fprintf('Percent of rejected draws: %.2f%%\n',...
        Prior.FractionBadDraws*100);
fprintf('log-prior correction: %.6f\n',Prior.LogTruncationCorrection);

p.PriorMean = mean(xd,2);
p.PriorSD = std(xd,0,2);
for jPrc=1:nPrc
    p.(PrcList{jPrc}) = prctile(xd,100*obj.ParamPercentiles(jPrc),2);
end
% Prior.Param = p;
obj.Param = p;

xdAux(:,BadDraws) = [];
pAux = obj.AuxParam;
pAux.PriorMean = mean(xdAux,2);
pAux.PriorSD = std(xdAux,0,2);
for jPrc=1:nPrc
    pAux.(PrcList{jPrc}) = prctile(xdAux,100*obj.ParamPercentiles(jPrc),2);
end
obj.AuxParam = pAux;
% Prior.AuxParam = pAux;

obj.Prior = Prior;

pList = {'Param','AuxParam'};
DispList = {'    Mean','PriorMean';
            '      SD','PriorSD';
            '      5%','PriorPrc050';
            '  Median','PriorPrc500';
            '     95%','PriorPrc950';
           }';
nc = size(DispList,2);
for jP=1:length(pList)
    Pj = pList{jP};
    psj = obj.(Pj);
    namelength = [cellfun('length',obj.(Pj).Names)];
    namelengthmax = max(namelength);
    fprintf(['\n%-',int2str(namelengthmax),'s'],'');
    for jc=1:nc
        fprintf('  %-8s',DispList{1,jc});
    end
    fprintf('\n');
    for jp=1:obj.(Pj).N
        fprintf(['%',int2str(namelengthmax),'s'],obj.(Pj).Names{jp});
        for jc=1:nc
            fprintf('  %8.4f',obj.(Pj).(DispList{2,jc})(jp));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Make Prior Report

fprintf('Making report: %s\n',ReportFileName);
fid = vcCreateTex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
fprintf(fid,'\\section{Parameters}\n');
str = [' & %.',int2str(obj.TablePrecision),'f'];
TableBreaks = obj.SetTableBreaks(obj.Param.N);
idxPar = 0;
nBreaks = length(TableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):TableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    if obj.TableMoveLeft
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
        fprintf(fid,'%s',obj.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Param.PriorDist{jr});
        fprintf(fid,str,Prior.UnconstrainedParam.PriorMode(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorMean(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorSD(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc050(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc500(jr));
        fprintf(fid,str,Prior.UnconstrainedParam.PriorPrc950(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Param.PriorMean(jr));
        fprintf(fid,str,obj.Param.PriorPrc050(jr));
        fprintf(fid,str,obj.Param.PriorPrc500(jr));
        fprintf(fid,str,obj.Param.PriorPrc950(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,obj.TableLines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\section{Auxiliary Parameters}\n');
str = [' & %.',int2str(obj.TablePrecision),'f'];
TableBreaks = obj.SetTableBreaks(obj.AuxParam.N);
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
        fprintf(fid,'%s',obj.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.AuxParam.PriorMean(jr));
        fprintf(fid,str,obj.AuxParam.PriorPrc050(jr));
        fprintf(fid,str,obj.AuxParam.PriorPrc500(jr));
        fprintf(fid,str,obj.AuxParam.PriorPrc950(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,obj.TableLines) && jr~=idxPar(end)
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
pdflatex(ReportFileName)

%% -------------------------------------------------------------------

%% Finish up
obj = obj.TrackTime(action,0);

end

function f = igamsolve(a,pmean,psd)
    f = 1./(a-1).*(pmean*gamma(a)./gamma(a-1/2)).^2-pmean^2-psd^2;
end
