function mcmcredux(obj,varargin)

% mcmcredux
% 
% Combine MCMC sample chains and save part of it. saved file contains a 
% structure named "draws" with fields:
%   draw.N (number of draws)
%   draw.Param (compact matrix with draws for Param)
%   draw.LPDF (row matrix with draws for posterior log-pdf)
%   draw.AuxParam (if op.Draws.AuxParam=1, matrix with draws for AuxParam)
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% Options
op.Draws.BurnIn = 0.25;
op.Draws.AuxParam = 0;
op.Draws.CombineChains = 1;
op.Draws.ExpandParam = 0;
op.NDraws = 10000;

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('Generating MCMC Draws Redux for Sample %.0f\n',obj.Post.MCMCStage)
ttName = sprintf('ReduxMCMC%.0f',obj.Post.MCMCStage);
obj.TimeTracker.start(ttName)

sample = obj.Post.MCMCSample;

nDrawsAvailable = floor(sample.NDraws*(1-op.Draws.BurnIn))*sample.NChains;
nDraws = min(op.NDraws,nDrawsAvailable);

%% load the mcmc draws
op.Draws.Thinning = max(1,...
    floor(sample.NDraws*sample.NChains*(1-op.Draws.BurnIn)/nDraws));
draws = obj.loadmcmcdraws(op.Draws);

%% Save MCMC draws Redux
fn = sprintf('%s_MCMC_%.0f_Redux',obj.Name,obj.Post.MCMCStage);
obj.Post.MCMCSample.FileNameRedux = fn;
save(fn,'-struct','draws')
fprintf('Saved MCMC draws redux to: %s.mat\n',fn)


%% Finish up
obj.TimeTracker.stop(ttName)

