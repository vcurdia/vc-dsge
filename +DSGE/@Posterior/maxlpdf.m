function maxlpdf(obj,varargin)

% maxlpdf
% 
% Search for the posterior mode by maximizing the log-pdf function
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 24, 2017
% Copyright (C) 2017 Vasco Curdia

%% Preamble
fprintf('\n*** Find posterior mode\n')
ttName = 'MaxLPDF';
obj.TimeElapsed.start(ttName)


%% settings
ReportFileName = sprintf('%s_Report_Posterior_Mode',obj.Model.Name);
ReportTitle = sprintf('%s\\\\Posterior Mode',obj.Model.Name);

np = obj.NEstimate;
pIdx = obj.EstimateIdx;

%% Options
op.NMax = 1;
op.ShowRobustness = 1;
op.KeepLogsMaxPost = 1;
op.KeepMatMaxPost = 0;
op.DrawAll=0;
op.PostArgIn={};
op.MinParams.verbose = 1;
op.MinParams.Guess.x0 = obj.Mode(pIdx);
MinParams.H0 = obj.Var(pIdx,pIdx);
op.GuessMaxDraws = 1000;
% op.GuessUsePriorDist = 0;
op.GuessPrcUsePriorDist = 0.5;
op.GuessMean = obj.Prior.Mean(pIdx);
op.GuessSD = obj.Prior.SD(pIdx);
op.GuessDF = 4;
op.Table = DSGE.Options.Table;


op = updateoptions(op,varargin);

%% postupdate checks and preparations
nx0 = size(op.MinParams.Guess.x0,2);
nMax = max(op.NMax,nx0);
nDrawPrior = floor((nMax-nx0)*op.GuessPrcUsePriorDist);
pNames = obj.Model.Param.Names(pIdx);
pNameLength = [cellfun('length',pNames)];
pNameLengthMax = max(pNameLength);

lpdfneg = @(x,varargin)(-obj.lpdf(x,varargin{:}));

jobOptions = cell(nMax,1);
for jm=1:nMax
    MinParamsj = op.MinParams;
    if jm<=nx0
        x0 = MinParamsj.Guess.x0(:,jm);
    else
        for jg=1:op.GuessMaxDraws
            if jm<nx0+nDrawPrior
                x0 = obj.Prior.Draw(1);
                x0 = x0(pIdx);
            else
                for jp=1:np
                    x0 = op.GuessMean;
                    for jpg=1:op.GuessMaxDraws*10
                        x0j = op.GuessMean(jp)+op.GuessSD(jp)*trnd(op.GuessDF);
                        if ismember(obj.Prior.Dist{jp},{'N'})
                            break
                        elseif ismember(obj.Prior.Dist{jp},{'B'}) ...
                                && x0j>0 && x0j<1
                            break
                        elseif ismember(obj.Prior.Dist{jp},...
                                        {'TN','G','IG1','IG2'}) && x0j>0
                            break
                        end
                    end
                    x0(jp,1) = x0j;
                end
            end
            f0 = lpdfneg(x0,PostArgIn{:});
            if f0<1e50
                break
            end
        end
        if f0<1e50
            error('Could not find acceptable guess for maximization %.0f.',jm)
        end
    end
    MinParamsj.Guess.x0 = x0;
    MinParamsj.MatFn = sprintf('MaxPost_%03.0f',jm);
    jobOptions{jm} = {lpdfneg,MinParamsj,sprintf('MaxPost_%03.0f.log',jm),...
                      op.PostArgIn};
end


%% Run minimizations
MaxPostOut = cell(1,nMax);
parfor jm=1:nMax
    MaxPostOut{jm} = maxlpdffcn(jobOptions{jm}{:});
end
MaxPostOut = [MaxPostOut{:}];
if ~op.KeepMatMaxPost
    for jm=1:nMax
        delete([jobOptions{jm}{2}.MatFn,'.mat'])
    end
end
if ~op.KeepLogsMaxPost
    for jm=1:nMax
        delete(jobOptions{jm}{3})
    end
end
nMax = length(MaxPostOut);

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

%% extract the best one
[ModeLPDF, idxMax] = min([MaxPostOut(:).f]);
obj.Mode(pIdx) = MaxPostOut(idxMax).x
obj.ModeLPDF = -MaxPostOut(idxMax).f
obj.Var(pIdx,pIdx) = MaxPostOut(idxMax).H;
obj.SD(pIdx) = diag(obj.Var).^(1/2);

%% display results on screen
pNames = obj.Model.Param.Names;
pNameLength = [cellfun('length',pNames)];
pNameLengthMax = max(pNameLength);
fprintf('\nResults from maximization of posterior:')
fprintf('\n=======================================\n')
DispList = {'','',pNames;
            'Prior','Dist',obj.Prior.Dist;
            '','  Mode',obj.Prior.Mode;
            '','  Mean',obj.Prior.Mean;
            '','   SD',obj.Prior.SD;
            '','   5%',obj.Prior.Prc05;
            '',' Median',obj.Prior.Median;
            '','   95%',obj.Prior.Prc95;
            'Posterior',' Mode ',obj.Mode;
            '','  SD',obj.SD;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(pNameLengthMax),'s'],DispList{1,jr});
    str2show = sprintf('%s  %-4s',str2show,DispList{2,jr});
    for jc=3:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for j=1:obj.Model.Param.N
    str2show = sprintf(['%',int2str(pNameLengthMax),'s'],DispList{1,3}{jp});
    str2show = sprintf('%s  %4s',str2show,DispList{2,3}{jp});
    for jc=3:nc
        str2show = sprintf('%s  %7.4f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\nposterior log-pdf at mode: %.6f\n\n',obj.ModeLPDF)

%% save minimization output
save([obj.Model.Name,'_MaxPostOut'],'MaxPostOut','idxMax')

%% create report
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);

% value of the posterior density
fprintf(fid,'\\begin{equation*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'posterior log density at mode: & %.4f\n',obj.ModeLPDF);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{equation*}\n');

fprintf(fid,'\\newpage \n');
fprintf(fid,'\\section{Parameters}\n');
np = obj.Model.Param.N;
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
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+7+1+2));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{7}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{2}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & Mode & Mean & SD & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mode & SD \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Prior.Dist{jr});
        fprintf(fid,str,obj.Prior.Mode(jr));
        fprintf(fid,str,obj.Prior.Mean(jr));
        fprintf(fid,str,obj.Prior.SD(jr));
        fprintf(fid,str,obj.Prior.Prc05(jr));
        fprintf(fid,str,obj.Prior.Median(jr));
        fprintf(fid,str,obj.Prior.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Mode(jr));
        fprintf(fid,str,obj.SD(jr));
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
np = obj.Model.AuxParam.N;
Mats = obj.Model.mats(obj.Mode,'SolveREE',0);
xAux = Mats.AuxParam;
str = [' & %.',int2str(op.Table.Precision),'f'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\section{Auxiliary Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+4+1+1));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior} & & Posterior\\\\[0.5ex]\n');
    fprintf(fid,'& Mean & 5\\%% & Median & 95\\%% & & Mode \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Sample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Sample.AuxParam.Prc95(jr));
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



%% Finish up
obj.TimeElapsed.stop(ttName)
end

function maxlpdffcn(fh,MinParams,logfn,PostArgIn)
    fid = fopen(logfn,'wt');
    PostArgIn = [{fid},PostArgIn];
    MinParams.varargin = PostArgIn;
    Out = vcrobustmin(fh,x0,MinParams,fid);
    if fid~=1,fclose(fid);end
end
