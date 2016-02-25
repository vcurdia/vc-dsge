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

s = CheckOptions(s);

if isfield(s.Options,'IRF')
    op = s.Options.IRF;
else
    op = struct; 
end

% Import default options from Sim
opList = fieldnames(s.Options.Sim);
for jOp=1:length(opList)
    Opj = opList{jOp};
    if ~isfield(op,Opj), op.(Opj) = s.Options.Sim.(opList{jOp}); end
end

% Check options
op.nPanels = length(op.Panels);
for jP=1:op.nPanels
    op.Panels(jP).nVar = length(op.Panels(jP).Var);
    if ~isfield(op.Panels(jP),'NamesPretty')
        op.Panels(jP).VarPretty = op.Panels(jP).Var;
    end
    if ~isfield(op.Panels(jP),'Scale') || isempty(op.Panels(jP).Scale)
        op.Panels(jP).Scale = ones(1,op.Panels(jP).nVar);
    end
end
op.nShocks2Show = length(op.Shocks2Show);

if ~isfield(op,'nSteps'), op.nSteps = 25; end
if ~isfield(op,'TickStep'), op.TickStep = 4; end

if ~isfield(op,'ShockSize'), op.ShockSize = ones(1,op.nShocks2Show); end

Fig = s.Options.Fig;
if isfield(op,'Fig')
    opList = fieldnames(op.Fig);
    for jOp=1:length(opList)
        Opj = opList{jOp};
        if ~isfield(Fig,Opj), Fig.(Opj) = op.Fig.(opList{jOp}); end
    end
end

if ~isfield(s,'PlotDir') || ~isfield(s.PlotDir,'IRF')
    s.PlotDir.IRF = 'Plots_IRF/';
end
if ~isdir(s.PlotDir.IRF), mkdir(s.PlotDir.IRF), end
s.FileName.PlotsIRF = sprintf('%s_IRF_%s',s.Spec,op.UseDist); 
op.ReportFilename = sprintf('%s_Report_IRF_%s',s.Spec,op.UseDist);
op.ReportTitle = sprintf('IRF Report:\\\\%s, %s',s.Spec,op.UseDist);

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
IRFObsVar = nan(s.n.ObsVar,op.nSteps,op.nShocks2Show,op.nDraws);
IRFStateVar = nan(s.n.StateVar,op.nSteps,op.nShocks2Show,op.nDraws);
IRFAuxVar = nan(s.n.AuxVar,op.nSteps,op.nShocks2Show,op.nDraws);
ShockIdx = zeros(op.nShocks2Show,1);
for j = 1:op.nShocks2Show
    ShockIdx(j) = find(ismember(s.ShockVar,op.Shocks2Show(j)));
end    
parfor jd=1:op.nDraws
    matj = fnMats(xd(:,jd));
    checkj = all(matj.REE.eu==1);
    if ~checkj
        IRFCheck(jd) = 0;
        continue
    end
    irf = zeros(s.n.StateVar,op.nShocks2Show,op.nSteps);
    irf = matj.REE.G2(:,ShockIdx);
    for t=2:op.nSteps
        irf(:,:,t) = matj.REE.G1*irf(:,:,t-1);
    end
    IRFStateVar(:,:,:,jd) = permute(irf,[1,3,2]);
    irfObs = zeros(s.n.ObsVar,op.nShocks2Show,op.nSteps);
    for t=1:op.nSteps
        irfObs(:,:,t) = matj.ObsEq.H*irf(:,:,t);
    end
    IRFObsVar(:,:,:,jd) = permute(irfObs,[1,3,2]);
    irfAux = zeros(s.n.AuxVar,op.nShocks2Show,op.nSteps);
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
for jS=1:op.nShocks2Show
    Sj = op.Shocks2Show{jS};
    for jP = 1:op.nPanels
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
        Figj = Fig;
        Figj.TitleList = Pj.VarPretty;
        Figj.FigShape = Pj.FigShape;
        OutFigj = vcFigure(PlotData,Figj);
        vcPrintPDF([s.PlotDir.IRF,s.FileName.PlotsIRF,...
             '_',Pj.Title,'_',Sj],Fig.KeepEPS,Fig.OpenPDF)
    end
end

%% Make report with IRF
fprintf('Making report: %s\n',op.ReportFilename);
fid = vcCreateTex(op.ReportFilename,op.ReportTitle);
for jS=1:op.nShocks2Show
    Sj = op.Shocks2Show{jS};
    fprintf(fid,'\\section{Shock: %s}\n',Sj);
    for jP = 1:op.nPanels
        Pj = op.Panels(jP).Title;
        fprintf(fid,'\\subsection{%s}\n',strrep(Pj,'_',': '));
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{IRF_%s_%s}\n',Pj,Sj);
        fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s_%s.pdf}\n',...
                s.PlotDir.IRF,s.FileName.PlotsIRF,Pj,Sj);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end
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
