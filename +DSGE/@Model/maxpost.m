function maxpost(obj,varargin)

% maxpost
% 
% Search for the posterior mode by maximizing the log-pdf function
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: March 24, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% Preamble
fprintf('\nSearching for Posterior Mode\n')
tmpFN = sprintf('%s-maxpost-tmp',obj.Name);
save(tmpFN)


%% settings
ReportFileName = sprintf('report-%s-param-postmode',obj.Name);
ReportTitle = sprintf('%s\\\\[30pt]Parameter Analysis\\\\Posterior Mode',...
                      obj.Name);

np = obj.Post.NEstimate;
pIdx = obj.Post.EstimateIdx;

%% Options
op.NMax = 10;
op.ShowRobustness = 1;
op.DrawAll = 0;
op.Min.verbose = 1;
op.Min.H0 = obj.Post.Var(pIdx,pIdx);
op.Guess = [obj.Post.Mode(pIdx),obj.Prior.Mean(pIdx)];
op.GuessMaxDraws = 1000;
% op.GuessUsePriorDist = 0;
op.GuessPrcUsePriorDist = 0.5;
op.GuessMean = obj.Prior.Mean(pIdx);
op.GuessSD = obj.Prior.SD(pIdx);
op.GuessDF = 4;
op.Table = DSGE.Options.Table;
op.KeepLogs = 1;
op.KeepMats = 1;

op = updateoptions(op,varargin{:});

%% postupdate checks and preparations
[npGuess,nx0] = size(op.Guess);
if npGuess==obj.Param.N
    op.Guess = op.Guess(pIdx,:);
end
nMax = max(op.NMax,nx0);
nDrawPrior = floor((nMax-nx0)*op.GuessPrcUsePriorDist);
pNames = obj.Param.Names(pIdx);
pNameLength = [cellfun('length',pNames)];
pNameLengthMax = max(pNameLength);
pDist = obj.Prior.Dist(pIdx);
lpdfneg = @(x,varargin)(-obj.postlpdf(x,varargin{:}));

x0 = zeros(np,nMax);
for jm=1:nMax
    if jm<=nx0
        x0(:,jm) = op.Guess(:,jm);
    else
        for jg=1:op.GuessMaxDraws
            if jm<nx0+nDrawPrior
                x0g = obj.priordraw(1);
                x0(:,jm) = x0g(pIdx);
            else
                for jp=1:np
                    for jpg=1:op.GuessMaxDraws*10
                        x0j = op.GuessMean(jp)+op.GuessSD(jp)*trnd(op.GuessDF);
                        if ismember(pDist{jp},{'N'})
                            break
                        elseif ismember(pDist{jp},{'B'}) && x0j>0 && x0j<1
                            break
                        elseif ismember(pDist{jp},{'TN','G','IG1','IG2'}) ...
                                && x0j>0
                            break
                        end
                    end
                    x0(jp,jm) = x0j;
                end
            end
            f0 = lpdfneg(x0(:,jm));
            if f0<1e50
                break
            end
        end
%         if f0>1e50
%             error('Could not find acceptable guess for maximization %.0f.',jm)
%         end
    end
end
save(tmpFN)

%% Run minimizations
MaxPostOut = cell(1,nMax);
x0j = op.GuessMean;
parfor jm=1:nMax
    fid = fopen(sprintf('%s-maxpost-%03.0f.log',obj.Name,jm),'wt');
    opj = op.Min;
    opj.MatFn = sprintf('%s-maxpost-%03.0f',obj.Name,jm);
    opj.LogFn = fid;
    opj.varargin = {struct('verbose',op.Min.verbose,'fid',fid)};
    MaxPostOut{jm} = robustmin(lpdfneg,x0(:,jm),opj);
    fclose(fid);
end
MaxPostOut = [MaxPostOut{:}];
nMax = length(MaxPostOut);
save(tmpFN)

%% extract the best one
[LPDFMode, idxMax] = min([MaxPostOut(:).f]);
obj.Post.Mode(pIdx) = MaxPostOut(idxMax).x;
obj.Post.ModeLPDF = -MaxPostOut(idxMax).f;
obj.Post.Var(pIdx,pIdx) = MaxPostOut(idxMax).H;
obj.Post.SD = diag(obj.Post.Var).^(1/2);

Mats = obj.mats(obj.Post.Mode,'SolveREE',0);
xAux = Mats.AuxParam;


%% save minimization output
save([obj.Name,'-maxpost-out'],'MaxPostOut','idxMax')

%% Show history evolution of robustness
if op.ShowRobustness
    for jm=1:nMax
        outj = MaxPostOut(jm).MinOutput;
        fprintf('\nRobustness analysis for minimization %3.0f',jm)
        fprintf('\n----------------------------------------\n')
        if isempty(outj)
            fprintf('%s\n%s\n\n',MaxPostOut(jm).rcMsg,MaxPostOut(jm).RrcMsg)
            continue
        end
        fprintf('Iteration %3.0f: function value: %15.8f',0,outj(1).fh)
%         fprintf('         %15s','')
        fprintf(' change: %15.8f',outj(1).fh-MaxPostOut(jm).f0)
        fprintf(' stopped at iteration %.0f\n',outj(1).itct)
        itbest.fh = outj(1).fh;
        itbest.idx = 1;
        for jr=2:length(outj)
            fprintf('Iteration %3.0f: function value: %15.8f',jr-1,outj(jr).fh)
            itchg = outj(jr).fh-itbest.fh;
            fprintf(' change: %15.8f',itchg)
            fprintf(' stopped at iteration %.0f\n',outj(jr).itct)
            if itchg<0
                itbest.fh = outj(jr).fh;
                itbest.idx = jr;
            end
        end
        fprintf('Best iteration: %.0f\n',itbest.idx-1)
        fprintf('Best iteration function value: %.8f\n',itbest.fh)
        fprintf('Robustness change over initial min: %.8f\n',...
                itbest.fh-outj(1).fh)
        fprintf('Robustness message: %s\n',MaxPostOut(jm).RrcMsg)
        fprintf('Best iteration message: %s\n\n',outj(itbest.idx).rcMsg)
    end
end

%% Show Starting Values
fprintf('\nGuess Values used (same order as in MaxPostOut):')
fprintf('\n================================================\n\n')
str2show = '    ';
for jp=1:np
    str2show = sprintf(['%s %',int2str(max(pNameLength(jp),9)),'s'],...
                       str2show,pNames{jp});
end
disp(str2show)
for jm=1:nMax
    str2show = sprintf('%4.0f',jm);
    for jp=1:np
        str2show = sprintf(['%s %',int2str(max(pNameLength(jp),9)),'.4f'],...
                           str2show,MaxPostOut(jm).x0(jp));
    end
    disp(str2show)
end
disp(' ')

%% show results for each maximization
fprintf('\nIndividual maximization results (ordered from best to worst):')
fprintf('\n=============================================================\n\n')
str2show = '     log-density';
for jp=1:np
    str2show = sprintf(['%s %',int2str(max(pNameLength(jp),9)),'s'],...
                       str2show,pNames{jp});
end
disp(str2show)
[SortPost,idxSortPost] = sort([MaxPostOut(:).f]);
for jm=1:nMax
    jShow = idxSortPost(jm);
    str2show = sprintf('%4.0f %11.4f',jShow,-MaxPostOut(jShow).f);
    for jp=1:np
        str2show = sprintf(['%s %',int2str(max(pNameLength(jp),9)),'.4f'],...
                           str2show,MaxPostOut(jShow).x(jp));
    end
    disp(str2show)
end
disp(' ')

%% show results on screen
pNames = obj.Param.Names;
pNameLength = [cellfun('length',pNames)];
nameLengthMax = max(pNameLength);
auxN = obj.AuxParam.N;
auxNames = obj.AuxParam.Names;
auxNameLength = [cellfun('length',auxNames)];
nameLengthMax = max(nameLengthMax,max(auxNameLength));

fprintf('\nResults from maximization of posterior:')
fprintf('\n=======================================\n')
DispList = {'','',pNames;
            'Prior','Dist',obj.Prior.Dist;
            '','   Mode',obj.Prior.Mode;
            '','   Mean',obj.Prior.Mean;
            '','     SD',obj.Prior.SD;
            '','     5%',obj.Prior.Prc05;
            '',' Median',obj.Prior.Median;
            '','    95%',obj.Prior.Prc95;
            'Posterior','   Mode',obj.Post.Mode;
            '','     SD',obj.Post.SD;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
    str2show = sprintf('%s  %-5s',str2show,DispList{2,jr});
    for jc=3:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:obj.Param.N
    str2show = sprintf(['%',int2str(nameLengthMax),'s'],DispList{1,3}{jp});
    str2show = sprintf('%s  %5s',str2show,DispList{2,3}{jp});
    for jc=3:nc
        str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\n')

DispList = {'','',auxNames;
            'Prior','   Mean',obj.Prior.Sample.AuxParam.Mean;
            '','     5%',obj.Prior.Sample.AuxParam.Prc05;
            '',' Median',obj.Prior.Sample.AuxParam.Median;
            '','    95%',obj.Prior.Sample.AuxParam.Prc95;
            'Posterior','   Mode',xAux;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
    for jc=2:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:auxN
    str2show = sprintf(['%',int2str(nameLengthMax),'s'],DispList{1,3}{jp});
    for jc=2:nc
        str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\n')

fprintf('\nposterior log-pdf at mode: %.6f\n\n',obj.Post.ModeLPDF)

%% create report
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);

fprintf(fid,'\\begin{equation*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'posterior log density at mode: & %.4f\n',obj.Post.ModeLPDF);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{equation*}\n');

fprintf(fid,'\\newpage \n');
np = obj.Param.N;
str = [' & $%.',int2str(op.Table.Precision),'f$'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\subsection{Parameters}\n');
    else
        fprintf(fid,'\\subsection{Parameters (%.0f/%.0f)}\n',...
                jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    if op.Table.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+7+1+2));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{7}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{2}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mode & SD \n');
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
        fprintf(fid,str,obj.Post.Mode(jr));
        fprintf(fid,str,obj.Post.SD(jr));
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

tableBreaks = settablebreaks(auxN,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\subsection{Auxiliary Parameters}\n');
    else
        fprintf(fid,'\\subsection{Auxiliary Parameters (%.0f/%.0f)}\n',...
                jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+4+1+1));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior} & & Posterior\\\\[0.5ex]\n');
    fprintf(fid,'& Mean & 5\\%% & Median & 95\\%% & & Mode \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,xAux(jr));
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

for jm=1:nMax
    if ~op.KeepLogs
        delete(sprintf('%s-maxpost-%03.0f.log',obj.Name,jm))
    end
    if ~op.KeepMats
        delete(sprintf('%s-maxpost-%03.0f.mat',obj.Name,jm))
    end
end


%% Finish up
delete([tmpFN,'.mat'])
end

