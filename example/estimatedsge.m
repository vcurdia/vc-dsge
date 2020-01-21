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

%% Load DSGE
cd(specname)
load(specname)
fprintf('Model: %s\n\n',model.Name)

%% MaxPost
vcdiary(specname,'maxpost')
tt.start('maxpost')
op.MaxPost.NMax = 10; %10
% op.MaxPost.Min.nit = 10; %comment this line to use default (1000)
model.maxpost(op.MaxPost)
save([specname,'%s-maxpost'])
model.sim('Dist','PostMode')
tt.stop('maxpost')
diary off

%% MCMC
op.MCMC.NDraws = 50000; %50000
% op.MCMC.NDrawsCalibrate = 100; %comment this line to use default (1000)
for s=1
    vcdiary(specname,sprintf('mcmc-%.0f',s))
    tt.start(sprintf('mcmc%.0f',s))
    model.Post.MCMCStage = s;
    model.mcmc(op.MCMC)
    save(sprintf('%s-mcmc-%.0f',specname,s))
    model.sim('Dist','PostDraws')
    save(specname)
    tt.stop(sprintf('mcmc%.0f',s))
    diary off
end


%% Finish up
tt.show
save(specname)
cd ..

