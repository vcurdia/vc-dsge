function mcmcredux(obj,varargin)

% mcmcredux
% 
% Combine MCMC sample chains and extract part of it
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.BurnIn = 0.5;
op.NDrawsRedux = 10000;

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Generating MCMC Draws Redux for Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('ReduxMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

sample = obj.MCMCSample;

nDrawsAvailable = floor(sample.NDraws*(1-op.BurnIn))*sample.NChains;
nDrawsRedux = min(op.NDrawsRedux,nDrawsAvailable);

%% load the mcmc draws
ThinningRedux = max(1,...
    floor(sample.NDraws*sample.NChains*(1-op.BurnIn)/nDrawsRedux));
idxDraws = (op.BurnIn*sample.NDraws+1):ThinningRedux:sample.NDraws;
for jChain=1:sample.NChains
    load(sample.FileNameDraws{jChain})
    xd(:,:,jChain) = xDraws(:,idxDraws);
    lpdfd(:,:,jChain) = lpdfDraws(:,idxDraws);
end
nDrawsUsed = size(xd,2)*sample.NChains;
xDraws = reshape(xd,obj.NEstimate,nDrawsUsed);
lpdfDraws = reshape(lpdfd,1,nDrawsUsed);
fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', ThinningRedux)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
obj.MCMCSample.NDrawsRedux = nDrawsUsed;

%% Save MCMC draws Redux
fn = sprintf('%s_MCMC_%.0f_Redux',obj.Model.Name,obj.MCMCStage);
obj.MCMCSample.FileNameRedux = fn;
save(fn,'xDraws','lpdfDraws')
fprintf('Saved MCMC draws redux to: %s.mat\n',fn)



%% Finish up
obj.TimeElapsed.stop(ttName)

