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
tt = TimeTracker;

%% Load DSGE
specName = 'MyDSGE';
specPath = specName;
basePath = cd(specPath);
load(specName)

%% Initiate parallel pool
% parpool(2)


%% MaxPost
op.MaxPost.NMax = 4;
op.MaxPost.Min.nit = 100;
op.MaxPost.Min.Ritmax = 4;
op.MaxPost.Min.Ritmin = 2;
post.maxlpdf(op.MaxPost)
save([specName,'_MaxPost'])
model.sim(post.Mode,op.Sim,'FNSuffix','_PostMode')


%% MCMC
op.MCMC.NDraws = 1000;
op.MCMC.NDrawsCalibrate = 200;
for s=1
    post.MCMCStage = s;
    post.mcmc(op.MCMC)
    save([specName,'_MCMC_',int2str(s)])
    save(specName)
    model.sim(post.draw(100),op.Sim,'FNSuffix','_PostDraws')
end


%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),tt.show

% exit


