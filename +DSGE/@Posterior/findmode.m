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

%% Options
op.NMAx = 1;
op.ShowRobustness = 1;
op.KeepLogsMaxPost = 1;
op.KeepMatMaxPost = 0;
op.DrawAll=0;
op.PrcGuessUsePriorDist=0.5;
op.PostArgIn={};
op.MinParams.verbose = 1;
op.MinParams.Guess.x0 = obj.Model.Param.Values(obj.Prior.EstimateIdx);

op = updateoptions(op,varargin);

%% postupdate checks
nx0 = size(op.MinParams.Guess.x0,2);
nMax = max(nMax,nx0);

HERE HERE HERE

JobOptions = cell(nMax,1);
for jm=1:nMax
    if exist('MinParams','var'), MinParamsj = MinParams; end
    if jm<=nx0
        MinParamsj.Guess.x0 = MinParamsj.Guess.x0(:,jm);
    else
        MinParamsj.Guess = rmfield(MinParamsj.Guess,'x0');
    end
    if jm<=nx0+floor((nMax-nx0)*PrcGuessUsePriorDist);
        MinParamsj.Guess.UsePriorDist = 1;
    end
    JobOptions{jm} = {FileName.Post,Params,jm,MinParamsj,...
                      sprintf('MaxPost_%03.0f.log',jm),PostArgIn};
end

%% other settings
ReportFileName = sprintf('%s_Report_Posterior_Mode',obj.Model.Name);
ReportTitle = sprintf('%s\\\\Posterior Mode',obj.Model.Name);

