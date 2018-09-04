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
op.BurnIn = 0.25;
op.Thinning = 1;
op.AuxParam = 0;
op.CombineChains = 1;
op.ExpandParam = 1;

op = updateoptions(op,varargin{:});

fprintf('Loading MCMC draws from sample %.0f\n',obj.MCMCStage)

sample = obj.MCMCSample;
draws.N = 0;

idxDraws = (op.BurnIn*sample.NDraws+1):op.Thinning:sample.NDraws;
for jChain=1:sample.NChains
    dc = load(sample.FileNameDraws{jChain});
    draws.Param(:,:,jChain) = dc.Param(:,idxDraws);
    draws.LPDF(:,:,jChain) = dc.LPDF(:,idxDraws);
end
nDraws = size(draws.LPDF,2);

if op.AuxParam
    fprintf('Generating AuxParam draws\n')
    dAux = zeros(obj.Model.AuxParam.N,nDraws,sample.NChains);
    fh = @(x)obj.Model.mats(obj.expandparam(x),'SolveREE',0);
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
    draws.Param = reshape(draws.Param,obj.NEstimate,draws.N);
    draws.LPDF = reshape(draws.LPDF,1,draws.N);
    if op.AuxParam
        draws.AuxParam = reshape(draws.AuxParam,obj.Model.AuxParam.N,draws.N);
    end
end

if op.ExpandParam
    draws.Param = obj.expandparam(draws.Param);
end

fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', op.Thinning)
fprintf('Total number of draws used: %.0f\n', draws.N)
fprintf('\n')
