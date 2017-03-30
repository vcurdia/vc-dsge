function mcmc(obj,varargin)

% mcmc
% 
% Generate MCMC sample
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 30, 2017
% Copyright (C) 2017 Vasco Curdia

%% Preamble
fprintf('\n*** Make MCMC sample, update %.0f\n',obj.MCMCUpdate)
ttName = sprintf('MCMCUpdate%.0f',obj.MCMCUpdate);
obj.TimeElapsed.start(ttName)


%% settings
% ReportFileName = sprintf('%s_Report_Posterior_Mode',obj.Model.Name);
% ReportTitle = sprintf('%s\\\\Posterior Mode',obj.Model.Name);

np = obj.NEstimate;
pIdx = obj.EstimateIdx;

fn = sprintf('%s_MCMC_Update_%.0f',obj.Model.Name,obj.MCMCUpdate); 

%% Options
op.KeepLogs = 1;
op.Chain.UseOneInitDrawAtMode = 0;

op = updateoptions(op,varargin{:});



%% Finish up
obj.TimeElapsed.stop(ttName)


