function obj = MakeIRF(obj,op)

% MakeIRF
% 
% Generates IRF for DSGE model
% 
% See also:
% DSGE, SetupMyDSGE
%
% .............................................................................
% 
% Created: November 17, 2016 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia

%-----------------------------------------------------------------------------%

%% Preamble

tt = TimeTracker;

fprintf('\n*** Making IRF\n')

%% Default Options
op.FigPanels = [];
op.Shocks2Show = obj.ShockVar.Names;
op.ShockSize = [];
op.Dist = 'Values';
op.nDraws = 1;
op.PlotDir = 'Plots_IRF/';
op.KeepPlots = 0;

%% Check options
obj = obj.CheckFigPanels;
nShocks2Show = length(op.Shocks2Show);
if isempty(op.ShockSize), op.ShockSize = ones(1,nShocks2Show); end

mkdir(op.PlotDir)
PlotFileName = sprintf('%s_IRF_%s',obj.Name,op.Dist); 
ReportFileName = sprintf('%s_Report_IRF_%s',obj.Name,op.Dist);
ReportTitle = sprintf('IRF Report:\\\\%s, %s',obj.Name,op.Dist);

%% -------------------------------------------------------------------

%% Prepare Draws
xd = obj.Param.GenDraws(op.Dist,op.NDraws);
nDraws = size(xd,2);

%% Generate IRF
fprintf('Generating IRFs...\n');
fnmats = @(x)obj.Mats(x,...
               'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
IRFCheck = ones(1,nDraws);
nSteps = obj.IRFNSteps;
nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;
IRF = nan(nStateVar+nObsVar+nAuxVar,nSteps,nShocks2Show,nDraws);
ShockIdx = zeros(nShocks2Show,1);
for j = 1:nShocks2Show
    ShockIdx(j) = find(ismember(obj.ShockVar.Names,obj.Shocks2Show(j)));
end   
parfor jd=1:nDraws
    matj = fnmats(xd(:,jd));
    checkj = all(matj.REE.eu==1);
    if ~checkj
        IRFCheck(jd) = 0;
        continue
    end
    irf = zeros(nStateVar,nShocks2Show,nSteps);
    irf(:,:,1) = matj.REE.G2(:,ShockIdx);
    for t=2:nSteps
        irf(:,:,t) = matj.REE.G1*irf(:,:,t-1);
    end
    IRFj = irf;
    if nObsVar>0
        irfObs = zeros(nObsVar,nShocks2Show,nSteps);
        for t=1:nSteps
            irfObs(:,:,t) = matj.ObsEq.H*irf(:,:,t);
        end
        IRFj = [IRFj;irfObs];
    end
    if nAuxVar>0
        irfAux = zeros(nAuxVar,nShocks2Show,nSteps);
        irfAux(:,:,1) = matj.AuxREE.G2(:,ShockIdx);
        for t=2:nSteps
            irfAux(:,:,t) = matj.AuxREE.G1*irf(:,:,t-1);
        end
        IRFj = [IRFj;irfAux];
    end
    IRF(:,:,:,jd) = permute(IRFj,[1,3,2]);
end
IRF(:,:,:,~IRFCheck) = [];
IRFCheck(~IRFCheck) = [];
nDrawsUsed = length(IRFCheck);

%% Plot IRFs
fprintf('Plotting IRFs...\n');
Fig = obj.Fig;
Fig.PlotBands = (nDrawsUsed>1);
Fig.XTick = 1:4:nSteps;
Fig.XTickLabel = 0:4:(nSteps-1);
VarNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
nPanels = length(obj.FigPanels);
for jP = 1:nPanels
    Pj = obj.FigPanels(jP);
    Figj = Fig;
    Figj.TitleList = Pj.PrettyNames;
    if isfield(Pj,'FigShape');
        Figj.FigShape = Pj.FigShape;
    end
    PlotData = nan(nDrawsUsed,nSteps,Pj.NVar,nShocks2Show);
    for jV=1:Pj.NVar
        Vj = Pj.Var{jV};
        [tf,idxV] = ismember(Vj,VarNames);
        if tf
            for jS=1:nShocks2Show
                PlotData(:,:,jV,jS) = obj.VarScale*Pj.Scale(jV)*...
                    obj.IRFShockSize(jS)*squeeze(IRF(idxV,:,jS,:))';
            end
        end
    end
    for jS=1:nShocks2Show
        if jS==1
            h = vcFigure(PlotData(:,:,:,jS),Figj);
        else
            h = vcFigureUpdate(h,PlotData(:,:,:,jS));
        end
        vcPrintPDF([obj.PlotDir.IRF,PlotFileName,...
             '_',Pj.Title,'_',obj.Shocks2Show{jS}],Fig.KeepEPS,Fig.OpenPDF)
    end
end

%% Make report with IRF
fprintf('Making report: %s\n',ReportFileName);
fid = vcCreateTex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
for jS=1:nShocks2Show
    Sj = obj.Shocks2Show{jS};
    fprintf(fid,'\\section{Shock: %s}\n',strrep(Sj,'_',''));
    for jP = 1:nPanels
        Pj = obj.FigPanels(jP).Title;
        fprintf(fid,'\\subsection{%s}\n',strrep(Pj,'_',': '));
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{IRF_%s_%s}\n',Pj,Sj);
        fprintf(fid,'\\includegraphics[scale=1]{%s%s_%s_%s.pdf}\n',...
                obj.PlotDir.IRF,PlotFileName,Pj,Sj);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)

%% -------------------------------------------------------------------

%% Finish up
close all
if ~op.KeepPlots, rmdir(op.PlotDir,'s'), end
tt.Show

end
%% -------------------------------------------------------------------
