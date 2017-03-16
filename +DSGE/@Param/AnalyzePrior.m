function AnalyzePrior(obj)

% AnalyzePrior
%
% Analyzes the priors
%
% See also:
% SetupMyDSGE, Param
%
% ...........................................................................
% 
% Created: November 10, 2016 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia


%% Preamble

fprintf('\n*** Analyzing DSGE Prior distribution\n')
TimeElapsed = tic;

%% Options

ReportFileName = sprintf('%s_Report_Prior',obj.Model.Name);
ReportTitle = sprintf('Prior Analysis:\\\\%s',obj.Model.Name);


%% Prepare variables
nPrc = length(obj.Percentiles);
PrcList = cell(nPrc,1);
for jPrc=1:nPrc
    PrcList{jPrc} = sprintf('Prc%03.0f',1000*obj.Percentiles(jPrc));
end

%% Analyze Parameters
np = obj.N;
obj.Prior.Mode = nan(np,1);
obj.Prior.Median = nan(np,1);
for jPrc=1:nPrc
    obj.Prior.(PrcList{jPrc}) = nan(np,1);
end
obj.Prior.DistParams = cell(np,1);
obj.Prior.LPdfCmd = cell(np,1);
obj.Prior.PdfCmd = cell(np,1);
obj.Prior.RndCmd = cell(np,1);
pOptions = optimoptions(@fsolve);
pOptions.Display = 'off';
for j=1:np

    if strcmp(obj.Prior.Dist{j},'C')
        obj.Prior.Median(j) = obj.Prior.Mean(j);
        obj.Prior.Mode(j) = obj.Prior.Mean(j);
        obj.Prior.SD(j) = 0;
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = obj.Prior.Mean(j);
        end
        obj.Prior.RndCmd{j} = @(n)repmat(obj.Prior.Mean(j),1,n);
    
    elseif strcmp(obj.Prior.Dist{j},'N')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        obj.Prior.Mode(j) = pmean;
        obj.Prior.Median(j) = pmean;
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = ...
                norminv(obj.Percentiles(jPrc),pmean,psd);
        end
        obj.Prior.DistParams{j} = [pmean,psd];
        obj.Prior.LPdfCmd{j} = sprintf('log(normpdf(%s,%.16f,%.16f))',...
                                    obj.Names{j},pmean,psd);
        obj.Prior.PdfCmd{j} = sprintf('normpdf(%%.16f,%.16f,%.16f)',pmean,psd);
        obj.Prior.RndCmd{j} = @(n)normrnd(pmean,psd,1,n);
    
    elseif strcmp(obj.Prior.Dist{j},'TN')
        % Assume x>=0
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = -pmean/psd;
        acdf = normcdf(a,0,1);
        aZ = 1-normcdf(a,0,1);
        alambda = normpdf(a,0,1)/aZ;
        adelta = alambda*(alambda-a);
        obj.Prior.Mean(j) = pmean + psd*alambda;
        obj.Prior.SD(j) = psd*(1-adelta)^(1/2);
        obj.Prior.Mode(j) = max(0,pmean);
        obj.Prior.Median(j) = norminv(0.5*aZ+acdf,pmean,psd);
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = norminv(...
                obj.Percentiles(jPrc)*aZ+acdf,pmean,psd);
        end
        obj.Prior.DistParams{j} = [pmean,psd];
        obj.Prior.LPdfCmd{j} = sprintf(...
            ['log((%1$s>=0)'...
             '*normpdf((%1$s-%2$.16f)/%3$.16f,0,1)/%3$.16f/%4$.16f)'],...
            obj.Names{j},pmean,psd,aZ);
        obj.Prior.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>=0)',...
             '*normpdf((%%1$.16f-%1$.16f)/%2$.16f,0,1)/%2$.16f/%3$.16f'],...
            pmean,psd,aZ);
        obj.Prior.RndCmd{j} = @(n)norminv(rand(1,n).*aZ+zcdf,pmean,psd);
    
    elseif strcmp(obj.Prior.Dist{j},'B')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = pmean*(pmean-pmean^2-psd^2)/psd^2;
        b = a*(1/pmean-1);
        obj.Prior.Mode(j) = min(max(0,(a-1)/(a+b-2)),1);
        obj.Prior.Median(j) = betainv(0.5,a,b);
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = betainv(obj.Percentiles(jPrc),a,b);
        end
        obj.Prior.DistParams{j} = [a,b];
        obj.Prior.LPdfCmd{j} = sprintf('log(betapdf(%s,%.16f,%.16f))',...
                                    obj.Names{j},a,b);
        obj.Prior.PdfCmd{j} = sprintf('betapdf(%%.16f,%.16f,%.16f)',a,b);
        obj.Prior.RndCmd{j} = @(n)betarnd(a,b,1,n);
        
    elseif strcmp(obj.Prior.Dist{j},'G')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = (pmean/psd)^2;
        b = pmean/a;
        if a>=1
            obj.Prior.Mode(j) = (a-1)*b;
        else
            obj.Prior.Mode(j) = NaN;
        end
        obj.Prior.Median(j) = gaminv(0.5,a,b);
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = gaminv(obj.Percentiles(jPrc),a,b);
        end
        obj.Prior.DistParams{j} = [a,b];
        obj.Prior.LPdfCmd{j} = sprintf('log(gampdf(%s,%.16f,%.16f))',...
                                    obj.Names{j},a,b);
        obj.Prior.PdfCmd{j} = sprintf('gampdf(%%.16f,%.16f,%.16f)',a,b);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b,1,n);
        
    elseif strcmp(obj.Prior.Dist{j},'IG1')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
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
        obj.Prior.Mode(j) = (1/b/(a+1/2))^(1/2);
        obj.Prior.Median(j) = gaminv(0.5,a,b)^(-1/2);
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = ...
                gaminv(1-obj.Percentiles(jPrc),a,b)^(-1/2);
        end
        obj.Prior.DistParams{j} = [a,b];
        obj.Prior.LPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-2),%2$.16f,%3$.16f)*2/%1$s^3))',...
            obj.Names{j},a,b);
        obj.Prior.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-2),%1$.16f,%2$.16f)*2/%%1$.16f^3)'],a,b);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b,1,n).^(-1/2);
    
    elseif strcmp(obj.Prior.Dist{j},'IG2')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        if psd==inf
            a = 2;
        else
            a = 2+pmean^2/psd^2;
        end
        b = 1/pmean/(a-1);
        obj.Prior.Mode(j) = 1/b/(a+1);
        obj.Prior.Median(j) = gaminv(0.5,a,b)^(-1);
        for jPrc=1:nPrc
            obj.Prior.(PrcList{jPrc})(j) = ...
                gaminv(1-obj.Percentiles(jPrc),a,b)^(-1);
        end
        obj.Prior.DistParams{j} = [a,b];
        obj.Prior.LPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',...
            p.Names{j},a,b);
        obj.Prior.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)'],a,b);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b).^(-1);
        
    end
end

%% display results on screen
fprintf('\nPrior:')
fprintf('\n------\n')
namelength = [cellfun('length',obj.Names)];
namelengthmax = max(namelength);
DispList = {'','Names';
            'Dist','Dist';
            '    Mode','Mode';
            '    Mean','Mean';
            '      SD','SD';
            '      5%','Prc050';
            '  Median','Median';
            '     95%','Prc950';
           }';
nc = size(DispList,2);
fprintf(['%-',int2str(namelengthmax),'s'],DispList{1,1});
fprintf('  %-4s',DispList{1,2});
for jc=3:nc
    fprintf('  %-8s',DispList{1,jc});
end
fprintf('\n');
for j=1:np
    fprintf(['%',int2str(namelengthmax),'s'],obj.(DispList{2,1}){j});
    fprintf('  %4s',obj.Prior.(DispList{2,2}){j});
    for jc=3:nc
        fprintf('  %8.4f',obj.Prior.(DispList{2,jc})(j));
    end
    fprintf('\n');
end
fprintf('\n');


%% Generate prior draws

fprintf('\nPrior Sample:')
fprintf('\n-------------\n\n')
xj = zeros(np,1);
BadDraws = false(1,obj.PriorNDraws);
xd = obj.DrawPrior(obj.PriorNDraws);
fh = @(x)obj.Model.Mats(x);
AuxNames = obj.Model.AuxParam.Names;
nAux = obj.Model.AuxParam.N;
xdAux = zeros(nAux,obj.PriorNDraws);
% Matsd = cell(obj.PriorNDraws);
parfor jd=1:obj.PriorNDraws
    Matsj = fh(xd(:,jd));
    BadDraws(jd) = ~Matsj.Status==1;
    for jp=1:nAux
        xdAux(jp,jd) = Matsj.AuxParam.(AuxNames{jp});
    end
%     Matsd{jd} = Matsj;
end
obj.Prior.Sample.NDraws = obj.PriorNDraws;
obj.Prior.Sample.NBadDraws = sum(BadDraws);
obj.Prior.Sample.FractionBadDraws = obj.Prior.Sample.NBadDraws/obj.PriorNDraws;
obj.Prior.LogTruncationCorrection = -log(1-obj.Prior.Sample.FractionBadDraws);
xd(:,BadDraws) = [];
obj.Prior.Sample.NDraws = size(xd,2);
fprintf('Number of accepted draws: %.0f\n',obj.Prior.Sample.NDraws);
fprintf('Percent of rejected draws: %.2f%%\n',...
        obj.Prior.Sample.FractionBadDraws*100);
fprintf('log-prior correction: %.6f\n',obj.Prior.LogTruncationCorrection);

obj.Prior.Sample.Param.Mean = mean(xd,2);
obj.Prior.Sample.Param.SD = std(xd,0,2);
obj.Prior.Sample.Param.Median = prctile(xd,50,2);
for jPrc=1:nPrc
    obj.Prior.Sample.Param.(PrcList{jPrc}) = ...
        prctile(xd,100*obj.Percentiles(jPrc),2);
end

xdAux(:,BadDraws) = [];
obj.Prior.Sample.AuxParam.Mean = mean(xdAux,2);
obj.Prior.Sample.AuxParam.SD = std(xdAux,0,2);
obj.Prior.Sample.AuxParam.Median = prctile(xdAux,50,2);
for jPrc=1:nPrc
    obj.Prior.Sample.AuxParam.(PrcList{jPrc}) = ...
        prctile(xdAux,100*obj.Percentiles(jPrc),2);
end

pList = {'Param','AuxParam'};
DispList = {'    Mean','Mean';
            '      SD','SD';
            '      5%','Prc050';
            '  Median','Median';
            '     95%','Prc950';
           }';
nc = size(DispList,2);
for jP=1:length(pList)
    Pj = pList{jP};
    psj = obj.Model.(Pj);
    namelength = [cellfun('length',obj.Model.(Pj).Names)];
    namelengthmax = max(namelength);
    fprintf(['\n%-',int2str(namelengthmax),'s'],'');
    for jc=1:nc
        fprintf('  %-8s',DispList{1,jc});
    end
    fprintf('\n');
    for jp=1:obj.Model.(Pj).N
        fprintf(['%',int2str(namelengthmax),'s'],obj.Model.(Pj).Names{jp});
        for jc=1:nc
            fprintf('  %8.4f',obj.Prior.Sample.(Pj).(DispList{2,jc})(jp));
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
TableBreaks = obj.SetTableBreaks(np);
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
    fprintf(fid,'& \\multicolumn{7}{c}{Prior Definition} ');
    fprintf(fid,'& & \\multicolumn{4}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mean & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Prior.Dist{jr});
        fprintf(fid,str,obj.Prior.Mode(jr));
        fprintf(fid,str,obj.Prior.Mean(jr));
        fprintf(fid,str,obj.Prior.SD(jr));
        fprintf(fid,str,obj.Prior.Prc050(jr));
        fprintf(fid,str,obj.Prior.Median(jr));
        fprintf(fid,str,obj.Prior.Prc950(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Prior.Sample.Param.Mean(jr));
        fprintf(fid,str,obj.Prior.Sample.Param.Prc050(jr));
        fprintf(fid,str,obj.Prior.Sample.Param.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.Param.Prc950(jr));
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
TableBreaks = obj.SetTableBreaks(nAux);
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
        fprintf(fid,'%s',obj.Model.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc050(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc950(jr));
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

%% Finish up
fprintf('\nAnalyzePrior: '), vctoc(TimeElapsed)

end

%% functions used

function f = igamsolve(a,pmean,psd)
    f = 1./(a-1).*(pmean*gamma(a)./gamma(a-1/2)).^2-pmean^2-psd^2;
end
