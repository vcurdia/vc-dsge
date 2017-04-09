function draws = loaddraws(obj,varargin)

% loaddraws
% 
% Load MCMC sample draws
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
op.Thinning = 1;
op.AuxParam = 0;

op = updateoptions(op,varargin{:});

fprintf('\nLoading MCMC draws from sample %.0f\n',obj.MCMCStage)

sample = obj.MCMCSample;
draws.N = 0;

idxDraws = (op.BurnIn*sample.NDraws+1):op.Thinning:obj.sample.NDraws;
for jChain=1:sample.NChains
    load(sample.FileNameDraws{jChain})
    xd(:,:,jChain) = xDraws(:,idxDraws);
    lpdfd(:,:,jChain) = lpdfDraws(:,idxDraws);
end
draws.N = size(xd,2)*sample.NChains;
draws.Param = obj.expandparam(reshape(xd,obj.NEstimate,draws.N));
draws.LPDF = reshape(lpdfd,1,draws.N);

if op.AuxParam
    fprintf('Generating AuxParam draws\n')
    xdAux = zeros(obj.Model.AuxParam.N,draws.N);
    fh = @(x)obj.Model.mats(x,'SolveREE',0);
    parfor jd=1:draws.N
        Matsj = fh(xd(:,jd));
        xdAux(:,jd) = Matsj.AuxParam;
    end
    draws.AuxParam = xdAux;
end

fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', op.Thinning)
fprintf('Total number of draws used: %.0f\n', draws.N)
