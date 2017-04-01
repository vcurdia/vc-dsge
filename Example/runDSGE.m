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
setpath
set(0,'defaultTextInterpreter','latex');
TimeElapsed = TimeTracker;

%% Settings
specName = 'MyDSGE';
specPath = specName;

%% Initiate parallel pool
% parpool(2)

%% Load DSGE
basePath = cd(specPath);
load(specName)
DSGE.linkobj(model,prior,post)

%% MaxPost
options.MaxPost.NMax = 4;
options.MaxPost.Min.nit = 2;
options.MaxPost.Min.Ritmax = 2;
options.MaxPost.Min.Ritmin = 1;
post.maxlpdf(options.MaxPost)
save(specName)
save([specName,'_MaxPost'])

%% MCMC
% post.MCMCNChains = 4;
% post.MCMCNDraws = 100000;
% post.MCMCBurnIn = 0.25;
% post.MCMCNThinning = 1;
% for post.MCMCUpdate=0
%     fprintf('\n*** MCMC Update %.0f *',post.MCMCUpdate)
%     post.JumpSize = 2.4;
%     post.calibratejumpsize
%     save(sprintf('%s_MCMC_Update%.0f_SSF',specName,post.MCMCUpdate))
%     post.mcmc
%     save(specName)
%     save(sprintf('%s_MCMC_Update%.0f',specName,post.MCMCUpdate))
%     delete(sprintf('%s_MCMC_Update%.0f_SSF.mat',specName,post.MCMCUpdate))
%     post.analyzemcmc
%     save(FileName.Output)
%     save(sprintf('%s_MCMC_Update%.0f',specName,post.MCMCUpdate))
% end



%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


