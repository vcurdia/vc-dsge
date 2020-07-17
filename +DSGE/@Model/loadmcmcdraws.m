function draws = loadmcmcdraws(obj,varargin)

% loadmcmcdraws
% 
% Load MCMC sample draws
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% Options
op.BurnIn = 0.25;
op.AuxParam = 0;
op.CombineChains = 1;
op.ExpandParam = 1;
op = updateoptions(op,varargin{:});

fprintf('Loading MCMC draws from sample %.0f\n',obj.Post.MCMCStage)

if ~isfield(obj.Post.MCMCSample,'NDrawsKeep')
    obj.Post.MCMCSample.NDrawsKeep = obj.Post.MCMCSample.NDraws; 
end

sample = obj.Post.MCMCSample;
draws.N = 0;

nDrawsKeep = min(sample.NDrawsKeep,sample.NDraws);
nThinning = sample.NDraws/sample.NDrawsKeep;
% idxDraws = (op.BurnIn*sample.NDraws+1):op.Thinning:sample.NDraws;
for jChain=1:sample.NChains
    dc = load(sample.FileNameDraws{jChain});
    nDrawsKeep = min(nDrawsKeep,size(draws.Param,2));
    idxDraws = (op.BurnIn*nDrawsKeep+1):nDrawsKeep;
    draws.Param(:,:,jChain) = dc.Param(:,idxDraws);
    draws.LPDF(:,:,jChain) = dc.LPDF(:,idxDraws);
end
nDraws = nDrawsKeep*nThinning;

if op.AuxParam
    fprintf('Generating AuxParam draws\n')
    dAux = zeros(obj.AuxParam.N,nDraws,sample.NChains);
    fh = @(x)obj.mats(obj.expandparam(x));
    for jC=1:sample.NChains
        parfor jd=1:nDraws
            Matsj = fh(draws.Param(:,jd,jC));
            dAux(:,jd,jC) = Matsj.AuxParam;
        end
    end
    draws.AuxParam = dAux;
end

draws.N = nDraws*sample.NChains;
if op.CombineChains
    draws.Param = reshape(draws.Param,obj.Post.NEstimate,draws.N);
    draws.LPDF = reshape(draws.LPDF,1,draws.N);
    if op.AuxParam
        draws.AuxParam = reshape(draws.AuxParam,obj.AuxParam.N,draws.N);
    end
end

if op.ExpandParam
    draws.Param = obj.expandparam(draws.Param);
end

fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', nThinning)
fprintf('Total number of draws used: %.0f\n', draws.N)
fprintf('\n')

end

