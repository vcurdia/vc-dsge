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
% ReportFileName = sprintf('%s_Report_Posterior_Mode',obj.Model.Name);
% ReportTitle = sprintf('%s\\\\Posterior Mode',obj.Model.Name);

sid = op.SampleID;
pIdx = obj.EstimateIdx;

fprintf('\n*** Make MCMC sample %.0f\n',sid)
ttName = sprintf('MCMCSample%.0f',sid);
obj.TimeElapsed.start(ttName)

obj.MCMCSample(sid).NChains = op.NChains;
obj.MCMCSample(sid).NDraws = op.Chain.NDraws;
if ~op.Chain.Augment
    op.Chain.JumpVar = op.JumpScale^2/obj.NEstimate*obj.Var(pIdx,pIdx);
    obj.MCMCSample(sid).JumpScale = op.JumpScale;
    obj.MCMCSample(sid).JumpVar = op.Chain.JumpVar;
    obj.MCMCSample(sid).NRejections = zeros(1,op.NChains);
else
    op.Chain.JumpVar = obj.MCMCSample(sid).JumpVar;
end
obj.MCMCSample(sid).FileNames = cell(op.NChains,1);
for jChain=1:op.NChains
    obj.MCMCSample(sid).FileNames{jChain} = sprintf(...
        '%s_MCMCSample_%.0f_Chain_%.0f',obj.Model.Name,sid,jChain);
end


%% Run MCMCFcn
[npx0,nx0] = size(op.x0);
if nx0>0 && npx0==obj.Model.Param.N, x0 = op.x0(obj.EstimateIdx,:);
nRejections = cell(op.NChains);
parfor jChain=1:op.NChains
    opj = op.Chain;
    opj.NRejections = obj.MCMCSample(sid).NRejections(jChain);
    opj.fn = obj.MCMCSample(sid).FileNames{jChain};
    if nx0>=jChain
        opj.x0 = x0(:,jChain);
    end
    nRejections{jChain} = obj.mcmcchain(opj);
end

    
%% show rejection rates
for jChain=1:op.NChains
    obj.MCMCSample(sid).NRejections(jChain) = nRejections{jChain};
    fprintf('Chain %.0f: JumpScaleFactor = %4.2f, Rejection rate = %5.1f%%\n',...
            jChain,op.Chain.JumpScaleFactor,...
            nRejections{jChain}/op.Chain.NDraws*100)
end
fprintf('\n')


%% Clean up
if ~op.KeepLogs
    for jChain=1:op.NChains
        delete(sprintf('%s.log',obj.MCMCSample(sid).FileNames{jChain},jChain));
    end
end

%% Finish up
obj.TimeElapsed.stop(ttName)


