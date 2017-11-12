% estimatedsge
%
% Example of how to estimate a DSGE using the VC-DSGE toolbox. It assumes that 
% the DSGE was previously setup by setupDSGE.
%
% See also:
% setupdsge, DSGE.Model, DSGE.Prior, DSGE.Posterior, DSGE.Data
%
% ...........................................................................
%
% Created: January 21, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia


%% Preamble
clear all
setpath
addpath('Mats/')
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
tt = TimeTracker;

%% Load DSGE
specName = 'MyDSGE';
specPath = specName;

cd(specPath)
load(specName)

%% Initiate parallel pool
% parpool(2)


%% MaxPost
op.MaxPost.NMax = 10;
op.MaxPost.Min.nit = 1000;
op.MaxPost.Min.Ritmax = 20;
op.MaxPost.Min.Ritmin = 5;
post.maxlpdf(op.MaxPost)
save([specName,'_MaxPost'])
model.sim(post.Mode,op.Sim,'FNSuffix','_PostMode')


%% MCMC
op.MCMC.NDraws = 50000;
op.MCMC.NDrawsCalibrate = 100;
for s=1
    post.MCMCStage = s;
    post.mcmc(op.MCMC)
    save([specName,'_MCMC_',int2str(s)])
    model.sim(post.draw(1000),op.Sim,'FNSuffix','_PostDraws')
    save(specName)
end


%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),tt.show

% exit


