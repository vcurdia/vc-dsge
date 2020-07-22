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
op.uselastdraw = 1;
op.UsePostDraw = 1;
op.NChains = 4;
op.JumpScale = 2.4;
op.JumpNewVarWeight = 1;
op.Augment = 0;
op.NDraws = 50000;
op.NDrawsKeep = 200000; 
op.CalibrateMCMC = [];
op.AnalyzePost = 1;
op.MCMCConvergence = 1;
op.MCMCRedux = 1;
op.Analysis = struct;
op = updateoptions(op,varargin{:});

tt = TimeTracker;

%% initial draw
if obj.Post.MCMCStage>1 && op.uselastdraw
    op.x0 = zeros(obj.Post.NEstimate,op.NChains);
    for jChain=1:obj.Post.MCMCSample.NChains
        draws = load(obj.Post.MCMCSample.FileNameDraws{jChain});
        op.x0(:,jChain) = draws.Param(:,end);
    end
    clear draws
end


%% MCMC calibration
if isempty(op.CalibrateMCMC)
    op.CalibrateMCMC = ~op.Augment;
end
if op.CalibrateMCMC
    tt.start('calibratemcmc')
    op.JumpScale = obj.calibratemcmc(op);
    tt.stop('calibratemcmc')
end

%% Preparations

if isempty(obj.Post.MCMCStage), obj.Post.MCMCStage = 1; end
fprintf('\nGenerating MCMC Sample %.0f\n',obj.Post.MCMCStage)
tmpFN = sprintf('%s-mcmc-%.0f-tmp',obj.Name,obj.Post.MCMCStage);
save(tmpFN)

pIdx = obj.Post.EstimateIdx;

[npx0,nx0] = size(op.x0);
x0 = cell(1,op.NChains);
if nx0>0
    if npx0==obj.Param.N
        op.x0 = op.x0(pIdx,:); 
    end
    for j=1:nx0
        x0{j} = op.x0(:,j);
    end
end
if nx0<op.NChains && op.UsePostDraw && obj.Post.MCMCStage>1
    x0d = obj.postdraw(op.NChains-nx0);
    for j=(nx0+1):op.NChains
        x0{j} = x0d(obj.Post.EstimateIdx,j-nx0);
    end
end

op.NDrawsKeep = min(op.NDrawsKeep,op.NDraws);


jumpVarRaw = obj.Post.Var(pIdx,pIdx)/obj.Post.NEstimate;
if op.JumpNewVarWeight<1 && obj.Post.MCMCStage>1
    jumpVarRaw = op.JumpNewVarWeight^2*jumpVarRaw + ...
        (1-op.JumpNewVarWeight)^2*obj.Post.MCMCSample.JumpVar/...
        obj.Post.MCMCSample.JumpScale^2;
end
obj.Post.MCMCSample.NChains = op.NChains;
obj.Post.MCMCSample.NDraws = op.NDraws;
obj.Post.MCMCSample.NDrawsKeep = op.NDrawsKeep;
obj.Post.MCMCSample.NThinning = op.NDraws/op.NDrawsKeep;
if ~op.Augment
    obj.Post.MCMCSample.JumpScale = op.JumpScale;
    obj.Post.MCMCSample.JumpVar = op.JumpScale^2*jumpVarRaw;
    obj.Post.MCMCSample.NRejections = zeros(1,op.NChains);
end
obj.Post.MCMCSample.FileNameDraws = cell(op.NChains,1);
for jChain=1:op.NChains
    obj.Post.MCMCSample.FileNameDraws{jChain} = sprintf(...
        '%s-mcmc-%.0f-chain-%.0f',obj.Name,obj.Post.MCMCStage,jChain);
end
obj.Post.MCMCSample.FileNameRedux = [];
save(tmpFN)

%% create MCMC chains
opChain.Augment = op.Augment;
opChain.NDraws = op.NDraws;
opChain.NDrawsKeep = op.NDrawsKeep;
opChain.JumpVar = obj.Post.MCMCSample.JumpVar;
nRejections = obj.Post.MCMCSample.NRejections;
tt.start('generatemcmc')
parfor jChain=1:op.NChains
    opj = opChain;
    opj.fn = obj.Post.MCMCSample.FileNameDraws{jChain}
    opj.x0 = x0{jChain};
    nRejections(jChain) = obj.mcmcchain(opj);
end
tt.stop('generatemcmc')
    
%% show rejection rates
obj.Post.MCMCSample.NRejections = nRejections;
for jChain=1:op.NChains
    fprintf('Chain %.0f: JumpScale = %4.2f, Rejection rate = %5.1f%%\n',...
            jChain,op.JumpScale,nRejections(jChain)/op.NDraws*100)
end
fprintf('\n')

%% Clean up
if ~op.KeepLogs
    for jChain=1:op.NChains
        delete(sprintf('%s.log',obj.Post.MCMCSample.FileNameDraws{jChain},jChain));
    end
end

%% Finish up MCMC
save(tmpFN)

%% Run MCMC analysis
if op.AnalyzePost
    tt.start('analyzepost')
    obj.analyzepost(op.Analysis)
    tt.stop('analyzepost')
    save(tmpFN)
end
if op.MCMCConvergence
    tt.start('mcmcconvergence')
    obj.mcmcconvergence(op.Analysis)
    tt.stop('mcmcconvergence')
    save(tmpFN)
end
if op.MCMCRedux, obj.mcmcredux, save(tmpFN), end
delete([tmpFN,'.mat'])

tt.showtimers



end

