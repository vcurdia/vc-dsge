function analyzeconvergence(obj,varargin)

% analyzeconvergence
% 
% Analyze convergence of MCMC sample
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 8, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.Draws.BurnIn = 0.5;
op.Draws.Thinning = 1;
op.Draws.AuxParam = 1;
op.Table = DSGE.Options.Table;
op.NBin = 50;
op.Fig = DSGE.Options.Figure;
op.Fig.Color = colorscheme;
op.Fig.FontSize = 6;
op.PlotDir = 'Plots_Convergence/';

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Analyzing convergence of MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('AnalyzeConvergenceMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
ReportFileName = sprintf('%s_Report_MCMC_%.0f_Convergence',obj.Model.Name,...
                         obj.MCMCStage);
ReportTitle = sprintf('%s\\\\MCMC Stage %.0f\\\\Convergence Analysis',...
                      obj.Model.Name,obj.MCMCStage);

sample = obj.MCMCSample;
np = obj.Model.Param.N;
pNames = obj.Model.Param.Names;
nAux = obj.Model.AuxParam.N;
auxNames = obj.Model.AuxParam.Names;

%% load the mcmc draws
draws = obj.loaddraws(op.Draws);





%% Finish up
obj.TimeElapsed.stop(ttName)

