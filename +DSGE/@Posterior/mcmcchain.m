function mcmcchain(obj,fn,nDraws,varargin)
function MCMCFcn(PostName,SaveName,Post,nDraws,jChain,options)

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
op.IsAugment = 0;
op.SaveList = {'xDraws','postDraws','nRejections','nDraws',...
               'jumpScaleFactor','jumpScale','jumpVar'};

op = updateoptions(op,varargin{:});

%% preparations
fid = fopen([fn,'.log'],'wt');
fpost = 

%% preliminary calculations
np = length(Post.Mode);
if ~exist('isAugment','var'), isAugment = 0; end


%% save output
save([fn,SaveList{:})

%% close printed output file
if fid~=1,fclose(fid);end

