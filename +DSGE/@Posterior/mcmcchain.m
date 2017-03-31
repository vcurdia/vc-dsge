function nRejections = mcmcchain(obj,varargin)

% mcmcchain
% 
% Generate MCMC chain
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 28, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.verbose = 1;
op.NDraws = 50000;
op.NIRS = 1000;
op.NBlocks = 10;
op.ExplosionScale = [.4,.5,.6,.7];
op.ExplosionProb = [.3,.6,.8,1];
op.NChangednIRS = 2;
op.InitDrawTStudent = 1;
op.InitDrawItMax = 1000;
op.InitDrawTol = 20;
op.InitDrawVarFactor = 1;
op.InitDrawDF = 6;
op.NRejections = 0;
op.fn = 'MCMC_Chain';

op = updateoptions(op,varargin{:});

%% preparations
fid = fopen([op.fn,'.log'],'wt');
np = obj.Model.Param.N;
nRejections = op.NRejections;



%% save output
save(op.fn,'xDraws','postDraws');

%% close printed output file
if fid~=1,fclose(fid);end

