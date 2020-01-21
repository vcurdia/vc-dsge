function analyzeprior(obj,varargin)

% analyzeprior
%
% Analyzes the prior distribution
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: November 10, 2016 by Vasco Curdia
% Copyright 2016-2018 by Vasco Curdia


%% Preamble
fprintf('\nAnalyzing Prior parameters\n')

%% Options
op.NDraws = 10000;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op.Table = DSGE.Options.Table;
op = updateoptions(op,varargin{:});

%% other settings
ReportFileName = sprintf('report-%s-param-prior',obj.Name);
ReportTitle = sprintf('%s\\\\[30pt]Parameter Analysis\\\\Prior',obj.Name);

%% useful variables
np = obj.Param.N;
pNames = obj.Param.Names;

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
    fprintf('  %4s',obj.Prior.(DispList{1,2}){j});
    for jc=2:nc
        fprintf('  %8.4f',obj.Prior.(DispList{jc,2})(j));
    end
    fprintf('\n');
end
fprintf('\n');


%% Generate prior draws
fprintf('\nPrior Sample:')
fprintf('\n-------------\n\n')
obj.makepriorsample(op)
AuxNames = obj.AuxParam.Names;
nAux = obj.AuxParam.N;

% Table
pList = {'Param','AuxParam'};
DispList = {'    Mean','Mean';
            '      SD','SD';
            '      5%','Prc05';
            '  Median','Median';
            '     95%','Prc95';
           };
nc = size(DispList,1);
namelengthmax = max([cellfun('length',{obj.Param.Names{:}, ...
                    obj.AuxParam.Names{:}})]);
for jP=1:length(pList)
    Pj = pList{jP};
    psj = obj.(Pj);
    fprintf(['\n%-',int2str(namelengthmax),'s'],'');
    for jc=1:nc
        fprintf('  %-8s',DispList{jc,1});
    end
    fprintf('\n');
    for jp=1:obj.(Pj).N
        fprintf(['%',int2str(namelengthmax),'s'],obj.(Pj).Names{jp});
        for jc=1:nc
            fprintf('  %8.4f',obj.Prior.Sample.(Pj).(DispList{jc,2})(jp));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Make Prior Report

fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
% fprintf(fid,'\\section{Parameters}\n');
str = [' & $%.',int2str(op.Table.Precision),'f$'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\section{Parameters}\n');
    else
        fprintf(fid,'\\section{Parameters (%.0f/%.0f)}\n',jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    if op.Table.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+7+1+3));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{7}{c}{Prior Definition} ');
    fprintf(fid,'& & \\multicolumn{3}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Prior.Dist{jr});
        fprintf(fid,str,obj.Prior.Mode(jr));
        fprintf(fid,str,obj.Prior.Mean(jr));
        fprintf(fid,str,obj.Prior.SD(jr));
        fprintf(fid,str,obj.Prior.Prc05(jr));
        fprintf(fid,str,obj.Prior.Median(jr));
        fprintf(fid,str,obj.Prior.Prc95(jr));
        fprintf(fid,' &');
         fprintf(fid,str,obj.Prior.Sample.Param.Prc05(jr));
        fprintf(fid,str,obj.Prior.Sample.Param.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.Param.Prc95(jr));
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

% fprintf(fid,'\\section{Auxiliary Parameters}\n');
str = [' & $%.',int2str(op.Table.Precision),'f$'];
tableBreaks = settablebreaks(nAux,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\section{Auxiliary Parameters}\n');
    else
        fprintf(fid,'\\section{Auxiliary Parameters (%.0f/%.0f)}\n',...
                jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+3));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{3}{c}{Prior Sample} \\\\[0.5ex]\n');
    fprintf(fid,'& 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc95(jr));
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

