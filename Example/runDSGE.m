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
tt = TimeTracker;

%% Load DSGE
specName = 'MyDSGE';
specPath = specName;
basePath = cd(specPath);
load(specName)
sim = DSGE.Sim(model,prior,post,data);

%% Initiate parallel pool
% parpool(2)

%% Settings
sim.TimeLabels = {'1990q1','1995q1','2000q1','2005q1'};


%% MaxPost
% options.MaxPost.NMax = 4;
% options.MaxPost.Min.nit = 100;
% options.MaxPost.Min.Ritmax = 10;
% options.MaxPost.Min.Ritmin = 5;
% post.maxlpdf(options.MaxPost)
% save([specName,'_MaxPost'])
% model.sim(data,post.Mode,options.Sim,'FNSuffix','_PostMode')

%% MCMC
options.MCMC.NDraws = 1000;
options.MCMC.NDrawsCalibrate = 200;
for s=1
    post.MCMCStage = s;
    post.mcmc(options.MCMC)
    save(specName)
    model.sim(data,post.draw(100),options.Sim,'FNSuffix','_PostDraws')
end


%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),tt.show

% exit


