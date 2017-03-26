function findmode(obj,varargin)

% findmode
% 
% Search for the posterior mode
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
ttName = 'FindMode';
obj.TimeElapsed.start(ttName)


%% Options
op.NMAx = 1;
op.ShowRobustness = 1;
op.KeepLogsMaxPost = 1;
op.KeepMatMaxPost = 0;
op.DrawAll=0;
op.PrcGuessUsePriorDist=0.5;
op.PostArgIn={};
op.MinParams.verbose = 1;
op.MinParams.Guess.x0 = obj.Mode;

op = updateoptions(op,varargin);

%% postupdate checks and preparations
nx0 = size(op.MinParams.Guess.x0,2);
nMax = max(nMax,nx0);
np = obj.Model.Param.N;
pNames = obj.Model.Param.Names;
pNameLength = [cellfun('length',pNames)];
pNameLengthMax = max(pNameLength);
pIdx = obj.Prior.EstimateIdx;

jobOptions = cell(nMax,1);
for jm=1:nMax
    MinParamsj = op.MinParams;
    if jm<=nx0
        MinParamsj.Guess.x0 = MinParamsj.Guess.x0(:,jm);
    else
        MinParamsj.Guess = rmfield(MinParamsj.Guess,'x0');
    end
    if jm<=nx0+floor((nMax-nx0)*op.PrcGuessUsePriorDist);
        MinParamsj.Guess.UsePriorDist = 1;
    end
    jobOptions{jm} = {jm,sprintf('MaxPost_%03.0f.log',jm),...
                      MinParamsj,op.PostArgIn};
end

%% other settings
ReportFileName = sprintf('%s_Report_Posterior_Mode',obj.Model.Name);
ReportTitle = sprintf('%s\\\\Posterior Mode',obj.Model.Name);

%% Run minimizations
MaxPostOut = cell(1,nMax);
parfor jm=1:nMax
    MaxPostOut{jm} = maxpost(jobOptions{jm}{:});
end
MaxPostOut = [MaxPostOut{:}];
if ~op.KeepMatMaxPost
    for jm=1:nMax
        delete(sprintf('MaxPost_%03.0f.mat',jm))
    end
end
if ~op.KeepLogsMaxPost
    for jm=1:nMax
        delete(jobOptions{jm}{2})
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
obj.Mode = MaxPostOut(idxMax).x
obj.ModeLPDF = -MaxPostOut(idxMax).f
obj.Var(pIdx,pIdx) = MaxPostOut(idxMax).H;
obj.SD(obj.Prior.EstimateIdx) = diag(obj.Var).^(1/2);

%% display results on screen
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
for j=1:np
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

%% create tex table
MakeTableMaxPost


%% Finish up
obj.TimeElapsed.stop(ttName)
