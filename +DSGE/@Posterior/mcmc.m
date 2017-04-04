function mcmc(obj,varargin)

% mcmc
% 
% Generate MCMC sample
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 30, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.KeepLogs = 1;
op.SampleID = length(obj.MCMCSample)+1;
op.x0 = [];
op.NChains = 4;
op.JumpScale = 2.4;
op.Chain.Augment = 0;
op.Chain.NDraws = 50000;

op = updateoptions(op,varargin{:});

%% Preparations

sid = op.SampleID;
pIdx = obj.EstimateIdx;

fprintf('\n*** Making MCMC Sample %.0f\n',sid)
ttName = sprintf('Sample%.0f',sid);
obj.TimeElapsed.start(ttName)

obj.Sample(sid).NChains = op.NChains;
obj.Sample(sid).NDraws = op.Chain.NDraws;
if ~op.Chain.Augment
    op.Chain.JumpVar = op.JumpScale^2/obj.NEstimate*obj.Var(pIdx,pIdx);
    obj.Sample(sid).JumpScale = op.JumpScale;
    obj.Sample(sid).JumpVar = op.Chain.JumpVar;
    obj.Sample(sid).NRejections = zeros(1,op.NChains);
else
    op.Chain.JumpVar = obj.Sample(sid).JumpVar;
end
obj.Sample(sid).FileNameDraws = cell(op.NChains,1);
for jChain=1:op.NChains
    obj.Sample(sid).FileNameDraws{jChain} = sprintf(...
        '%s_MCMC_Sample_%.0f_Chain_%.0f',obj.Model.Name,sid,jChain);
end
obj.Sample(sid).NDrawsRedux = [];
obj.Sample(sid).FileNameRedux = [];


%% create MCMC chains
[npx0,nx0] = size(op.x0);
x0 = cell(1,4);
if nx0>0
    if npx0==obj.Model.Param.N
        op.x0 = op.x0(obj.EstimateIdx,:); 
    end
    for j=1:nx0
        x0{j} = op.x0(:,j);
    end
end
nRejections = obj.Sample(sid).NRejections;
parfor jChain=1:op.NChains
    opj = op.Chain;
    opj.NRejections = nRejections(jChain);
    opj.fn = obj.Sample(sid).FileNameDraws{jChain}
    opj.x0 = x0{jChain};
    nRejections(jChain) = obj.mcmcchain(opj);
end

    
%% show rejection rates
for jChain=1:op.NChains
    obj.Sample(sid).NRejections(jChain) = nRejections(jChain);
    fprintf('Chain %.0f: JumpScale = %4.2f, Rejection rate = %5.1f%%\n',...
            jChain,op.JumpScale,nRejections(jChain)/op.Chain.NDraws*100)
end
fprintf('\n')


%% Clean up
if ~op.KeepLogs
    for jChain=1:op.NChains
        delete(sprintf('%s.log',obj.Sample(sid).FileNameDraws{jChain},...
                       jChain));
    end
end

%% Finish up
obj.TimeElapsed.stop(ttName)


