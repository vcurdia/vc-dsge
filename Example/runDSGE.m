% runMyDSGE
%
% This file gives an example of how to use the vcDSGE package
%
% See also:
% setupDSGE
%
% ...........................................................................
%
% Created: January 21, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia


%% Preamble
clear all
% startup
TimeElapsed = TimeTracker;

%% Load DSGE
specName = 'MyDSGE';
specPath = specName;
basePath = cd(specPath);
load(specName)
DSGE.linkobj(model,prior,post)

%% Initiate parallel pool
% parpool(2)

%% Settings
options.States.Tick.Labels = {'1990q1','1995q1','2000q1','2005q1'};


%% MaxPost
% options.MaxPost.NMax = 4;
% options.MaxPost.Min.nit = 100;
% options.MaxPost.Min.Ritmax = 10;
% options.MaxPost.Min.Ritmin = 5;
% post.maxlpdf(options.MaxPost)
% save([specName,'_MaxPost'])
% model.irf(post.Mode,'FileNameSuffix','_Mode')
% model.vd(post.Mode,'FileNameSuffix','_Mode')

%% MCMC
options.MCMC.NDraws = 1000;
options.MCMC.NDrawsCalibrate = 200;
for s=1
    post.MCMCStage = s;
    fnSuffix = sprintf('_MCMC_%.0f',s);
    fn = [specName,fnSuffix];
    options.MCMC.JumpScale = 2.4; %reset to default
    options.MCMC.JumpScale = post.calibratemcmc(options.MCMC);
    save([fn,'_CalibrateJump'])
    post.mcmc(options.MCMC)
    save(fn)
    post.analyzeparam
    post.mcmcredux
    xd = post.draw(100);
    model.irf(xd,'FileNameSuffix',fnSuffix)
    model.vd(xd,'FileNameSuffix',fnSuffix)
    model.states(data,xd,options.States,'FileNameSuffix',fnSuffix)
    save([fn,'_Analysis'])
    save(specName)
end


%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


