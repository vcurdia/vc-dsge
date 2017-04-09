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
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.Draws.BurnIn = 0.5;
op.Draws.AuxParam = 0;
op.NDraws = 10000;

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Generating MCMC Draws Redux for Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('ReduxMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

sample = obj.MCMCSample;

nDrawsAvailable = floor(sample.NDraws*(1-op.BurnIn))*sample.NChains;
nDraws = min(op.NDraws,nDrawsAvailable);

%% load the mcmc draws
op.Draws.Thinning = max(1,...
    floor(sample.NDraws*sample.NChains*(1-op.BurnIn)/nDraws));

draws = obj.loaddraws(op.Draws);
draws.Param = draws.Param(obj.EstimateIdx,:);

%% Save MCMC draws Redux
fn = sprintf('%s_MCMC_%.0f_Redux',obj.Model.Name,obj.MCMCStage);
obj.MCMCSample.FileNameRedux = fn;
save(fn,'draws')
fprintf('Saved MCMC draws redux to: %s.mat\n',fn)


%% Finish up
obj.TimeElapsed.stop(ttName)

