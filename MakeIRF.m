function s = MakeIRF(s)

% MakeIRF
% 
% Generates IRF for DSGE model
% 
% See also:
% MakeIRFFcn
%
% .............................................................................
% 
% Created: February 9, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%-----------------------------------------------------------------------------%

%% Preamble

Action = 'IRF';

fprintf('\n*** Making IRF\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% Options

if isfield(s.Options,'IRF')
    op = s.Options.IRF;
else
    op = struct; 
end

if ~isfield(op,'UseDist'), op.UseDist = 'PriorPrc500'; end
if ~isfield(op,'nDraws')
    if ismember(op.UseDist,{'PriorDraws','PostDraws'})
        op.nDraws = 1000;
    else
        op.nDraws = 1;
    end
end
if ~isfield(op,'Bands2Show'), op.Bands2Show = [50,70,90]; end
if ~isfield(op,'nPanelMax'), op.nPanelMax = 6; end
if ~isfield(op,'Vars2Show')
    list = {'ObsVar','StateVar','AuxVar'};
    jP = 0;
    for jlist=1:length(list)
        listj = list{jlist};
        for j=(1:ceil(s.n.(listj)/op.nPanelMax))
            jP = jP+1;
            op.Vars2Show(jP).Title = sprintf('%s%02.0f',listj,j);
            op.Vars2Show(jP).Names = s.(listj)(...
                (j-1)*op.nPanelMax+1:min(j*op.nPanelMax,s.n.(listj)));
        end
    end
end
if ~isfield(op,'Shocks2Show'), op.Shocks2Show = s.ShockVar; end

if ~isfield(op,'nSteps'), op.nSteps = 25; end
if ~isfield(op,'TickStep'), op.TickStep = 4; end

if ~isfield(op,'yPerSlack'), op.yPerSlack = 0.05; end
if ~isfield(op,'yMaxSlack'), op.yMaxSlack = []; end
if ~isfield(op,'yMinScale'), op.yMinScale = 0; end

if ~isfield(op,'ShowFig'), op.ShowFig = 0; end

if ~isfield(s,'PlotDir') || ~isfield(s.PlotDir,'IRF')
    s.PlotDir.IRF = 'Plots_IRF/';
end
if ~isdir(s.PlotDir.IRF), mkdir(s.PlotDir.IRF), end
if ~isfield(s.FileName,'PlotsIRF'), 
    s.FileName.PlotsIRF = sprintf('%sIRF',s.Spec); 
end

%% Check options
nPanels = length(op.Vars2Show);
for jP=1:nPanels
    if ~isfield(op.Vars2Show(jP),'NamesPretty')
        op.Vars2Show(jP).NamesPretty = op.Vars2Show(jP).Names;
    end
    if ~isfield(op.Vars2Show(jP),'Scale')
        op.Vars2Show(jP).Scale = 1;
    end
end
nShocks2Show = length(op.Shocks2Show);
if ~isfield(op,'ShockSize'), op.ShockSize = ones(1,nShocks2Show); end

% Save options
s.Options.IRF = op;

%% -------------------------------------------------------------------

%% Prepare Draws
if strcmp(op.UseDist,'PriorDraws')
    xd = feval(s.FileName.PriorDraw,op.nDraws);
elseif strcmp(op.UseDist,'PostDraws')
    load(FileName.MCMCDrawsRedux,'xd')
    xd = xd(:,1:op.nDraws);
else
    if ~isfield(s.Param,op.UseDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
        return
    end
    xd = s.Param.(op.UseDist);
end

%% Generate IRF
fprintf('Generating IRFs...\n');
fnMats = @(x)feval(s.FileName.Mats,x,...
               'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
IRFCheck = ones(1,op.nDraws);
IRFObsVar = nan(s.n.ObsVar,op.nSteps,nShocks2Show,op.nDraws);
IRFStateVar = nan(s.n.StateVar,op.nSteps,nShocks2Show,op.nDraws);
IRFAuxVar = nan(s.n.AuxVar,op.nSteps,nShocks2Show,op.nDraws);
ShockIdx = zeros(nShocks2Show,1);
for j = 1:nShocks2Show
    ShockIdx(j) = find(ismember(s.ShockVar,op.Shocks2Show(j)));
end    
parfor jd=1:op.nDraws
    matj = fnMats(xd(:,jd));
    checkj = all(matj.REE.eu==1);
    if ~checkj
        IRFCheck(jd) = 0;
        continue
    end
    irf = zeros(s.n.StateVar,nShocks2Show,op.nSteps);
    irf = matj.REE.G2(:,ShockIdx);
    for t=2:op.nSteps
        irf(:,:,t) = matj.REE.G1*irf(:,:,t-1);
    end
    IRFStateVar(:,:,:,jd) = permute(irf,[1,3,2]);
    irfObs = zeros(s.n.ObsVar,nShocks2Show,op.nSteps);
    for t=1:op.nSteps
        irfObs(:,:,t) = matj.ObsEq.H*irf(:,:,t);
    end
    IRFObsVar(:,:,:,jd) = permute(irfObs,[1,3,2]);
    irfAux = zeros(s.n.AuxVar,nShocks2Show,op.nSteps);
    irfAux(:,:,1) = matj.AuxREE.G2(:,ShockIdx);
    for t=2:op.nSteps
        irfAux(:,:,t) = matj.AuxREE.G1*irf(:,:,t-1);
    end
    IRFAuxVar(:,:,:,jd) = permute(irfAux,[1,3,2]);
end
keyboard
IRFObsVar(:,:,:,~IRFCheck) = [];
IRFStateVar(:,:,:,~IRFCheck) = [];
IRFAuxVar(:,:,:,~IRFCheck) = [];
IRFCheck(~IRFCheck) = [];
nDrawsUsed = length(IRFCheck);

keyboard


%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
