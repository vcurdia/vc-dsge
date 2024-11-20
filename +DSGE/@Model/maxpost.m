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
ReportFileName = sprintf('%s-param-postmode',obj.Name);
ReportTitle = sprintf('%s\\\\[30pt]Parameter Analysis\\\\Posterior Mode',...
                      obj.Name);

np = obj.Post.NEstimate;
pIdx = obj.Post.EstimateIdx;

%% Options
op.NMax = 20;
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

rng(0)
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
    rng(jm)
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

%% save minimization output
save([obj.Name,'-maxpost-out'],'MaxPostOut')

%% extract the best one
obj.maxpostchoosebest


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

for jm=1:nMax
    if ~op.KeepLogs
        delete(sprintf('%s-maxpost-%03.0f.log',obj.Name,jm))
    end
    if ~op.KeepMats
        delete(sprintf('%s-maxpost-%03.0f.mat',obj.Name,jm))
    end
end


%% make maxpost report
obj.maxpostreport(op)


%% Finish up
delete([tmpFN,'.mat'])

end

