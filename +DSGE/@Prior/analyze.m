function analyze(obj,varargin)

% analyze
%
% Analyzes a the prior distribution. Not
%
% See also:
% setupMyDSGE, DSGE.Prior
%
% ...........................................................................
% 
% Created: November 10, 2016 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia


%% Preamble
fprintf('\n*** Analyzing DSGE Prior\n')
ttName = 'Analyze';
obj.TimeElapsed.start(ttName)

%% Options
op.NDraws = 1000;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op.Table = DSGE.Options.Table;

op = updateoptions(op,varargin);

%% other settings
ReportFileName = sprintf('%s_Report_Prior',obj.Model.Name);
ReportTitle = sprintf('%s\\\\Prior Analysis',obj.Model.Name);

%% useful variables
np = obj.Model.Param.N;
pNames = obj.Model.Param.Names;

%% display results on screen
fprintf('\nPrior:')
fprintf('\n------\n')
namelength = [cellfun('length',pNames)];
namelengthmax = max(namelength);
DispList = {'Dist','Dist';
            '    Mode','Mode';
            '    Mean','Mean';
            '      SD','SD';
            '      5%','Prc05';
            '  Median','Median';
            '     95%','Prc95';
           };
nc = size(DispList,1);
fprintf(['%-',int2str(namelengthmax),'s'],'');
fprintf('  %-4s','');
for jc=2:nc
    fprintf('  %-8s',DispList{jc,1});
end
fprintf('\n');
for j=1:np
    fprintf(['%',int2str(namelengthmax),'s'],pNames{j});
    fprintf('  %4s',obj.(DispList{1,2}){j});
    for jc=2:nc
        fprintf('  %8.4f',obj.(DispList{jc,2})(j));
    end
    fprintf('\n');
end
fprintf('\n');


%% Generate prior draws

fprintf('\nPrior Sample:')
fprintf('\n-------------\n\n')
xj = zeros(np,1);
BadDraws = false(1,op.NDraws);
xd = obj.draw(op.NDraws);
fh = @(x)obj.Model.mats(x);
AuxNames = obj.Model.AuxParam.Names;
nAux = obj.Model.AuxParam.N;
xdAux = zeros(nAux,op.NDraws);
% Matsd = cell(obj.PriorNDraws);
parfor jd=1:op.NDraws
    Matsj = fh(xd(:,jd));
    BadDraws(jd) = ~Matsj.Status==1;
    xdAux(:,jd) = Matsj.AuxParam;
%     Matsd{jd} = Matsj;
end
obj.Sample.NDraws = op.NDraws;
obj.Sample.NBadDraws = sum(BadDraws);
obj.Sample.FractionBadDraws = obj.Sample.NBadDraws/op.NDraws;
obj.LPDFCorrection = -log(1-obj.Sample.FractionBadDraws);
xd(:,BadDraws) = [];
xdAux(:,BadDraws) = [];
obj.Sample.NDrawsUsed = size(xd,2);
fprintf('Number of accepted draws: %.0f\n',obj.Sample.NDrawsUsed);
fprintf('Percent of rejected draws: %.2f%%\n',...
        obj.Sample.FractionBadDraws*100);
fprintf('log-prior correction: %.6f\n',obj.LPDFCorrection);

obj.Sample.Param = sumstats(xd,op.Percentiles);
obj.Sample.AuxParam = sumstats(xdAux,op.Percentiles);

pList = {'Param','AuxParam'};
DispList = {'    Mean','Mean';
            '      SD','SD';
            '      5%','Prc05';
            '  Median','Median';
            '     95%','Prc95';
           };
nc = size(DispList,1);
for jP=1:length(pList)
    Pj = pList{jP};
    psj = obj.Model.(Pj);
    namelength = [cellfun('length',obj.Model.(Pj).Names)];
    namelengthmax = max(namelength);
    fprintf(['\n%-',int2str(namelengthmax),'s'],'');
    for jc=1:nc
        fprintf('  %-8s',DispList{jc,1});
    end
    fprintf('\n');
    for jp=1:obj.Model.(Pj).N
        fprintf(['%',int2str(namelengthmax),'s'],obj.Model.(Pj).Names{jp});
        for jc=1:nc
            fprintf('  %8.4f',obj.Sample.(Pj).(DispList{jc,2})(jp));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Make Prior Report

fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
fprintf(fid,'\\section{Parameters}\n');
str = [' & %.',int2str(op.Table.Precision),'f'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    if op.Table.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+7+1+5));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{7}{c}{Prior Definition} ');
    fprintf(fid,'& & \\multicolumn{4}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mean & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Dist{jr});
        fprintf(fid,str,obj.Mode(jr));
        fprintf(fid,str,obj.Mean(jr));
        fprintf(fid,str,obj.SD(jr));
        fprintf(fid,str,obj.Prc05(jr));
        fprintf(fid,str,obj.Median(jr));
        fprintf(fid,str,obj.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Sample.Param.Mean(jr));
        fprintf(fid,str,obj.Sample.Param.Prc05(jr));
        fprintf(fid,str,obj.Sample.Param.Median(jr));
        fprintf(fid,str,obj.Sample.Param.Prc95(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\section{Auxiliary Parameters}\n');
str = [' & %.',int2str(op.Table.Precision),'f'];
tableBreaks = settablebreaks(nAux,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Auxiliary Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+4));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Mean & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Sample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Prc95(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
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
obj.TimeElapsed.stop(ttName)

end

%% functions used

function f = igamsolve(a,pmean,psd)
    f = 1./(a-1).*(pmean*gamma(a)./gamma(a-1/2)).^2-pmean^2-psd^2;
end
