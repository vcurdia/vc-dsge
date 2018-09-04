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
op.NChains = 4;
op.JumpScale = 2.4;
op.Augment = 0;
op.NDraws = 50000;
op.CalibrateMCMC = [];
op.AnalyzeParam = 1;
op.AnalyzeConvergence = 1;
op.MCMCRedux = 1;


op = updateoptions(op,varargin{:});

%% MCMC calibration
if isempty(op.CalibrateMCMC)
    op.CalibrateMCMC = ~op.Augment;
end
if op.CalibrateMCMC
    op.JumpScale = obj.calibratemcmc(op);
    save(sprintf('_%s_MCMC_%.0f_CalibrateJump',obj.Model.Name,obj.MCMCStage))
end

%% Preparations

if isempty(obj.MCMCStage), obj.MCMCStage = 1; end
fprintf('\n*** Generating MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('MCMC%.0f',obj.MCMCStage);
obj.TimeTracker.start(ttName)

pIdx = obj.EstimateIdx;

obj.MCMCSample.NChains = op.NChains;
obj.MCMCSample.NDraws = op.NDraws;
if ~op.Augment
    obj.MCMCSample.JumpScale = op.JumpScale;
    obj.MCMCSample.JumpVar = op.JumpScale^2/obj.NEstimate*obj.Var(pIdx,pIdx);
    obj.MCMCSample.NRejections = zeros(1,op.NChains);
end
obj.MCMCSample.FileNameDraws = cell(op.NChains,1);
for jChain=1:op.NChains
    obj.MCMCSample.FileNameDraws{jChain} = sprintf(...
        '%s_MCMC_%.0f_Chain_%.0f',obj.Model.Name,obj.MCMCStage,jChain);
end
obj.MCMCSample.FileNameRedux = [];


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
opChain.Augment = op.Augment;
opChain.NDraws = op.NDraws;
opChain.JumpVar = obj.MCMCSample.JumpVar;
nRejections = obj.MCMCSample.NRejections;
parfor jChain=1:op.NChains
    opj = opChain;
    opj.fn = obj.MCMCSample.FileNameDraws{jChain}
    opj.x0 = x0{jChain};
    nRejections(jChain) = obj.mcmcchain(opj);
end
    
%% show rejection rates
obj.MCMCSample.NRejections = nRejections;
for jChain=1:op.NChains
    fprintf('Chain %.0f: JumpScale = %4.2f, Rejection rate = %5.1f%%\n',...
            jChain,op.JumpScale,nRejections(jChain)/op.NDraws*100)
end
fprintf('\n')

%% save workspace
save(sprintf('_%s_MCMC_%.0f',obj.Model.Name,obj.MCMCStage))

%% Clean up
if ~op.KeepLogs
    for jChain=1:op.NChains
        delete(sprintf('%s.log',obj.MCMCSample.FileNameDraws{jChain},jChain));
    end
end

%% Finish up MCMC
obj.TimeTracker.stop(ttName)

%% Run MCMC analysis
if op.AnalyzeParam, obj.analyzeparam, end
if op.AnalyzeConvergence, obj.analyzeconvergence, end
if op.MCMCRedux, obj.mcmcredux, end



