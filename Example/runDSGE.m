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

%% Initiate parallel pool
% parpool(2)

%% Settings
op.Sim.Data = data;
op.Sim.Tick.Labels = {'1990q1','1995q1','2000q1','2005q1'};


%% MaxPost
% op.MaxPost.NMax = 4;
% op.MaxPost.Min.nit = 100;
% op.MaxPost.Min.Ritmax = 10;
% op.MaxPost.Min.Ritmin = 5;
% post.maxlpdf(op.MaxPost)
% save([specName,'_MaxPost'])
% model.sim(post.Mode,op.Sim,'FNSuffix','_PostMode')

%% MCMC
op.MCMC.NDraws = 1000;
op.MCMC.NDrawsCalibrate = 200;
for s=1
    post.MCMCStage = s;
    post.mcmc(op.MCMC)
    save(specName)
    model.sim(post.draw(100),op.Sim,'FNSuffix','_PostDraws')
end


%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),tt.show

% exit


