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
options.MaxPost.Min.nit = 100;
options.MaxPost.Min.Ritmax = 10;
options.MaxPost.Min.Ritmin = 5;
post.maxlpdf(options.MaxPost)
save([specName,'_MaxPost'])

%% MCMC
options.MCMC.Chain.NDraws = 1000;
options.MCMC.NDrawsCalibrateJump = 100;
for sid=1
    options.MCMC.SampleID = sid;
    fn = sprintf('%s_MCMC_Sample_%.0f',specName,sid);
    options.MCMC.JumpScale = 2.4; %reset to default
    options.MCMC.JumpScale = post.calibratemcmcjump(options.MCMC);
    save([fn,'_CalibrateJump'])
    post.mcmc(options.MCMC)
    save(fn)
%     post.analyzemcmc
%     save([fn,'_Analysis'])
end



%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


