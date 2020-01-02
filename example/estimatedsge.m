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
specname = 'mydsge';

% Initiate parallel pool
%parpool('local',20)

%% Load DSGE
cd(specname)
load(specname)
fprintf('Model: %s\n\n',model.Name)

% %% MaxPost
% op.MaxPost.NMax = 10; %10
% op.MaxPost.Min.nit = 1; %1000;
% model.maxpost(op.MaxPost)
% save([specname,'-maxpost'])
% model.sim('Dist','PostMode')


%% MCMC
op.MCMC.NDraws = 100; %50000
op.MCMC.NDrawsCalibrate = 100; %1000
for s=1
    model.Post.MCMCStage = s;
    model.mcmc(op.MCMC)
    save([specname,'-mcmc-',int2str(s)])
    model.sim('Dist','PostDraws')
    save(specname)
end


%% Finish up
% delete(gcp)
fprintf('\n'),tt.show,fprintf('\n')
save(specname)
cd ..

% exit


