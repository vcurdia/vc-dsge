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
op.x0 = [];
op.UsePostDraw = 1;
op.NChains = 4;
op.JumpScale = 2.4;
op.Augment = 0;
op.NDraws = 50000;
op.CalibrateMCMC = [];
op.AnalyzePost = 1;
op.MCMCConvergence = 1;
op.MCMCRedux = 1;


op = updateoptions(op,varargin{:});

%% MCMC calibration
if isempty(op.CalibrateMCMC)
    op.CalibrateMCMC = ~op.Augment;
end
if op.CalibrateMCMC
    op.JumpScale = obj.calibratemcmc(op);
end

%% Preparations

if isempty(obj.Post.MCMCStage), obj.Post.MCMCStage = 1; end
fprintf('\n*** Generating MCMC Sample %.0f\n',obj.Post.MCMCStage)
ttName = sprintf('MCMC%.0f',obj.Post.MCMCStage);
obj.TimeTracker.start(ttName)

pIdx = obj.Post.EstimateIdx;

obj.Post.MCMCSample.NChains = op.NChains;
obj.Post.MCMCSample.NDraws = op.NDraws;
if ~op.Augment
    obj.Post.MCMCSample.JumpScale = op.JumpScale;
    obj.Post.MCMCSample.JumpVar = op.JumpScale^2/obj.Post.NEstimate*obj.Post.Var(pIdx,pIdx);
    obj.Post.MCMCSample.NRejections = zeros(1,op.NChains);
end
obj.Post.MCMCSample.FileNameDraws = cell(op.NChains,1);
for jChain=1:op.NChains
    obj.Post.MCMCSample.FileNameDraws{jChain} = sprintf(...
        '%s_MCMC_%.0f_Chain_%.0f',obj.Name,obj.Post.MCMCStage,jChain);
end
obj.Post.MCMCSample.FileNameRedux = [];


%% create MCMC chains
[npx0,nx0] = size(op.x0);
x0 = cell(1,op.NChains);
if nx0>0
    if npx0==obj.Param.N
        op.x0 = op.x0(obj.Post.EstimateIdx,:); 
    end
    for j=1:nx0
        x0{j} = op.x0(:,j);
    end
end
if op.UsePostDraw && obj.Post.MCMCStage>1
    x0d = obj.postdraw(op.NChains-nx0);
    for j=(nx0+1):op.NChains
        x0{j} = x0d(obj.Post.EstimateIdx,j-nx0);
    end
end

opChain.Augment = op.Augment;
opChain.NDraws = op.NDraws;
opChain.JumpVar = obj.Post.MCMCSample.JumpVar;
nRejections = obj.Post.MCMCSample.NRejections;
parfor jChain=1:op.NChains
    opj = opChain;
    opj.fn = obj.Post.MCMCSample.FileNameDraws{jChain}
    opj.x0 = x0{jChain};
    nRejections(jChain) = obj.mcmcchain(opj);
end
    
%% show rejection rates
obj.Post.MCMCSample.NRejections = nRejections;
for jChain=1:op.NChains
    fprintf('Chain %.0f: JumpScale = %4.2f, Rejection rate = %5.1f%%\n',...
            jChain,op.JumpScale,nRejections(jChain)/op.NDraws*100)
end
fprintf('\n')

%% save workspace
save(sprintf('%s_MCMC_%.0f',obj.Name,obj.Post.MCMCStage))

%% Clean up
if ~op.KeepLogs
    for jChain=1:op.NChains
        delete(sprintf('%s.log',obj.Post.MCMCSample.FileNameDraws{jChain},jChain));
    end
end

%% Finish up MCMC
obj.TimeTracker.stop(ttName)

%% Run MCMC analysis
if op.AnalyzePost, obj.analyzepost, end
if op.MCMCConvergence, obj.mcmcconvergence, end
if op.MCMCRedux, obj.mcmcredux, end



