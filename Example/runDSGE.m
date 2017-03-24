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

%% Settings
specName = 'MyDSGE';
specPath = specName;

%% Initiate parallel pool
% parpool(2)

%% Load DSGE
basePath = cd(specPath);
load(specName)
DSGE.linkobj(model,prior,post)

%% Finish up
% delete(gcp)
save(specName)
cd(basePath)
fprintf('\n'),TimeElapsed.show

% exit


