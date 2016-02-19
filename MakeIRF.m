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
if ~isfield(op,'Panels')
    list = {'ObsVar','StateVar','AuxVar'};
    jP = 0;
    for jlist=1:length(list)
        listj = list{jlist};
        for j=(1:ceil(s.n.(listj)/op.nPanelMax))
            jP = jP+1;
            op.Panels(jP).Title = sprintf('%s%02.0f',listj,j);
            op.Panels(jP).Var = s.(listj)(...
                (j-1)*op.nPanelMax+1:min(j*op.nPanelMax,s.n.(listj)));
        end
    end
end
if ~isfield(op,'Shocks2Show'), op.Shocks2Show = s.ShockVar; end
if ~isfield(op,'VarScale'), op.VarScale = 1; end

if ~isfield(op,'nSteps'), op.nSteps = 25; end
if ~isfield(op,'TickStep'), op.TickStep = 4; end

if ~isfield(s.Options,'Fig'), s.Options.Fig = struct; end
Fig = s.Options.Fig;
if ~isfield(Fig,'Visible'), Fig.Visible = 'off'; end
if ~isfield(Fig,'YSlack'), Fig.YSlack = 0.05; end
if ~isfield(Fig,'YMinScale'), Fig.YMinScale = 0; end
if ~isfield(Fig,'KeepEPS'), Fig.KeepEPS = 0; end
if ~isfield(Fig,'OpenPDF'), Fig.OpenPDF = 0; end
s.Options.Fig = Fig;

if ~isfield(s,'PlotDir') || ~isfield(s.PlotDir,'IRF')
    s.PlotDir.IRF = 'Plots_IRF/';
end
if ~isdir(s.PlotDir.IRF), mkdir(s.PlotDir.IRF), end
if ~isfield(s.FileName,'PlotsIRF'), 
    s.FileName.PlotsIRF = sprintf('%sIRF_%s',s.Spec,op.UseDist); 
end

if ~isfield(op,'ReportFilename')
    op.ReportFilename = sprintf('%sReport_IRF_%s',s.Spec,op.UseDist);
end
if ~isfield(op,'ReportTitle')
    op.ReportTitle = sprintf('IRF Report:\\\\%s, %s',s.Spec,op.UseDist);
end

%% Check options
nPanels = length(op.Panels);
for jP=1:nPanels
    op.Panels(jP).nVar = length(op.Panels(jP).Var);
    if ~isfield(op.Panels(jP),'NamesPretty')
        op.Panels(jP).VarPretty = op.Panels(jP).Var;
    end
    if ~isfield(op.Panels(jP),'Scale') || isempty(op.Panels(jP).Scale)
        op.Panels(jP).Scale = ones(1,op.Panels(jP).nVar);
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
IRFObsVar(:,:,:,~IRFCheck) = [];
IRFStateVar(:,:,:,~IRFCheck) = [];
IRFAuxVar(:,:,:,~IRFCheck) = [];
IRFCheck(~IRFCheck) = [];
nDrawsUsed = length(IRFCheck);

%% Plot IRFs
fprintf('Plotting IRFs...\n');
Fig.PlotBands = (nDrawsUsed>1);
Fig.XTick = 1:4:op.nSteps;
Fig.XTickLabel = 0:4:(op.nSteps-1);
for jS=1:nShocks2Show
    Sj = op.Shocks2Show{jS};
    for jP = 1:nPanels
        Pj = op.Panels(jP);
        PlotData = nan(nDrawsUsed,op.nSteps,Pj.nVar);
        for jV=1:Pj.nVar
            Vj = Pj.Var{jV};
            [tf,idxV] = ismember(Vj,s.ObsVar);
            if tf
               PlotData(:,:,jV) = Pj.Scale(jV)*op.ShockSize(jS)*...
                   squeeze(IRFObsVar(idxV,:,jS,:))';
            end
            if ~tf
                [tf,idxV] = ismember(Vj,s.StateVar);
                if tf
                    PlotData(:,:,jV) = Pj.Scale(jV)*op.ShockSize(jS)*...
                        squeeze(IRFStateVar(idxV,:,jS,:))';
                end
            end
            if ~tf
                [tf,idxV] = ismember(Vj,s.AuxVar);
                if tf
                    PlotData(:,:,jV) = Pj.Scale(jV)*op.ShockSize(jS)*...
                        squeeze(IRFAuxVar(idxV,:,jS,:))';
                end
            end
        end
        PlotData = op.VarScale*PlotData;
        Fig.TitleList = Pj.VarPretty;
        OutFigj = vcFigure(PlotData,Fig);
        vcPrintPDF(...
            [s.PlotDir.IRF,s.FileName.PlotsIRF,s.FileName.PlotsIRF,...
             '_',Pj.Title,'_',Sj],Fig.KeepEPS,Fig.OpenPDF)
    end
end

%% Make report with IRF
fprintf('Making report...\n');
fid = vcCreateTex(op.ReportFilename,op.ReportTitle);
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(op.ReportFilename)


%% -------------------------------------------------------------------

if strcmp(Fig.Visible,'off')
    close all
end

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);
fprintf('\n%s %s\n\n',Action,vctoc([],s.TimeElapsed.(Action)))

%% -------------------------------------------------------------------
