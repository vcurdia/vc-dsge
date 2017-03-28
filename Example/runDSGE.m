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
opMaxPost.NMAx = 1;
opMaxPost.MinParams.nit = 2; %1000
opMaxPost.MinParams.Ritmax = 1;%30;
opMaxPost.MinParams.Ritmin = 1;%10;
post.maxlpdf(opMaxPost)



%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


