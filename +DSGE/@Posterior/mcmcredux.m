function mcmcredux(obj,varargin)

% mcmcredux
% 
% Combine MCMC sample and extract part of it
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.SampleID = length(obj.MCMCSample);
op.BurnIn = 0.5;
op.NDrawsRedux = 10000;

op = updateoptions(op,varargin{:});

%% Preparations

sid = op.SampleID;
sample = obj.MCMCSample(sid);

fprintf('\n*** Generating MCMC Draws Redux for Sample %.0f\n',sid)
ttName = sprintf('MCMCSample%.0fRedux',sid);
obj.TimeElapsed.start(ttName)

nDrawsAvailable = floor(sample.NDraws*(1-op.BurnIn))*sample.NChains;
nDrawsRedux = min(op.NDrawsRedux,nDrawsAvailable);

%% load the mcmc draws
nThinningRedux = max(1,...
    floor(sample.NDraws*sample.NChains*(1-op.BurnIn)/nDrawsRedux));
idxDraws = (op.BurnIn*sample.NDraws+1):nThinningRedux:sample.NDraws;
for jChain=1:sample.NChains
    load(sample.FileNameDraws{jChain})
    xd(:,:,jChain) = xDraws(:,idxDraws);
    lpdfd(:,:,jChain) = lpdfDraws(:,idxDraws);
end
clear xDraws lpdfDraws
nDrawsUsed = size(xd,2)*sample.NChains;
xDraws = reshape(xd,obj.NEstimate,nDrawsUsed);
lpdfDraws = reshape(lpdfd,nDrawsUsed,1);
fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', nThinningRedux)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)
obj.MCMCSample(sid).NDrawsRedux = nDrawsUsed;

%% Save MCMC draws Redux
obj.MCMCSample(sid).FileNameRedux = sprintf(...
        '%s_MCMC_Sample_%.0f_Redux',obj.Model.Name,sid);
save(obj.MCMCSample(sid).FileNameRedux,'xDraws','lpdfDraws')
fprintf('Saved MCMC draws redux to: %s.mat\n',obj.MCMCSample(sid).FileNameRedux)



%% Finish up
obj.TimeElapsed.stop(ttName)

