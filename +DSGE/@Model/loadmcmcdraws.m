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
    op.silent = 0;
    op = updateoptions(op,varargin{:});

    if ~op.silent
        fprintf('Loading MCMC draws from sample %.0f\n',obj.Post.MCMCStage)
    end

    if ~isfield(obj.Post.MCMCSample,'NDrawsKeep')
        obj.Post.MCMCSample.NDrawsKeep = obj.Post.MCMCSample.NDraws; 
    end
    if ~isfield(obj.Post.MCMCSample,'NThinning')
        obj.Post.MCMCSample.NThinning = ...
            obj.Post.MCMCSample.NDraws/obj.Post.MCMCSample.NDrawsKeep; 
    end
    sample = obj.Post.MCMCSample;

    draws.N = 0;

    nDrawsKeep = sample.NDrawsKeep;
% idxDraws = (op.BurnIn*sample.NDraws+1):op.Thinning:sample.NDraws;
    for jChain=1:sample.NChains
        dc = load(sample.FileNameDraws{jChain});
        nDrawsKeep = min(nDrawsKeep,size(dc.Param,2));
        idxDraws = ceil((op.BurnIn*nDrawsKeep+1):nDrawsKeep);
        draws.Param(:,:,jChain) = dc.Param(:,idxDraws);
        draws.LPDF(:,:,jChain) = dc.LPDF(:,idxDraws);
    end
    nDraws = length(idxDraws);

    if op.AuxParam
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

    draws.N = nDraws;

    if ~op.CombineChains
        draws.NChains = sample.NChains;
        draws.NTotal = nDraws*sample.NChains;
    else
        draws.N = nDraws*sample.NChains;
        draws.Param = reshape(draws.Param,obj.Post.NEstimate,draws.N);
        draws.LPDF = reshape(draws.LPDF,1,draws.N);
        if op.AuxParam
            draws.AuxParam = reshape(draws.AuxParam,obj.AuxParam.N,draws.N);
        end
    end

    if op.ExpandParam
        draws.Param = obj.expandparam(draws.Param);
    end

    if ~op.silent
        fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
        fprintf('Number of chains: %.0f\n', sample.NChains)
        fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
        fprintf('Thinning used: %.1f\n', sample.NThinning)
        if op.CombineChains, fprintf('Chains were combined.\n'), end
        fprintf('Total number of draws used: %.0f\n', draws.N)
        fprintf('\n')
    end

end

