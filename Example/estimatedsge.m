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

% Spec
specName = 'LM';

% keep log
diary(sprintf('%s_%.0f-%02.0f-%02.0f-%02.0f%02.0f%02.0f_Estimation.log', ...
              specName,clock))
diary on

% Initiate parallel pool
%parpool('local',20)

%% Load DSGE
cd(specName)
load(specName)

%% MaxPost
op.MaxPost.NMax = 10; %10
op.MaxPost.Min.nit = 1000; %1000;
op.MaxPost.Min.Ritmax = 20; %20;
op.MaxPost.Min.Ritmin = 5;
model.maxpost(op.MaxPost)
save([specName,'_MaxPost'])
model.sim('Dist','PostMode')


%% MCMC
op.MCMC.NDraws = 50000;
op.MCMC.NDrawsCalibrate = 1000;
for s=1
    model.Post.MCMCStage = s;
    model.mcmc(op.MCMC)
    save([specName,'_MCMC_',int2str(s)])
    model.sim('Dist','PostDraws')
    save(specName)
end


%% Finish up
% delete(gcp)
save(specName)
tt.show
diary off
cd ..

% exit


