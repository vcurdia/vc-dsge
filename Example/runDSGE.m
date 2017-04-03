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
% options.MaxPost.NMax = 4;
% options.MaxPost.Min.nit = 2;
% options.MaxPost.Min.Ritmax = 2;
% options.MaxPost.Min.Ritmin = 1;
% post.maxlpdf(options.MaxPost)
% save(specName)
% save([specName,'_MaxPost'])

%% MCMC
options.MCMC.Chain.NDraws = 10;
for sid=1
    options.MCMC.SampleID = sid;
    options.MCMC.JumpScale = 2.4;
%     post.calibratejumpscale
%     save(sprintf('%s_MCMC_Sample_%.0f_SSF',specName,sid))
    post.mcmc(options.MCMC)
%     save(specName)
%     save(sprintf('%s_MCMC_Sample_%.0f',specName,sid))
%     delete(sprintf('%s_MCMC_Sample_%.0f_SSF.mat',specName,sid))
%     post.analyzemcmc
%     save(FileName.Output)
    save(sprintf('%s_MCMC_Sample_%.0f',specName,sid))
end



%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


