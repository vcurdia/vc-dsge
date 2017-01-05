function obj = MakeIRF(obj)

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
% Copyright 2016 by Vasco Curdia

%-----------------------------------------------------------------------------%

%% Preamble

action = 'IRF';
obj = obj.TrackTime(action,1);

fprintf('\n*** Making IRF\n')

%% Options

% Check options
obj = obj.CheckFigPanels;
nShocks2Show = length(obj.Shocks2Show);

if isempty(obj.IRFShockSize), obj.IRFShockSize = ones(1,nShocks2Show); end

if ~isfield(obj.PlotDir,'IRF')
    obj.PlotDir.IRF = 'Plots_IRF/';
end
mkdir(obj.PlotDir.IRF)
PlotFileName = sprintf('%s_IRF_%s',obj.Name,obj.SimDist); 
ReportFileName = sprintf('%s_Report_IRF_%s',obj.Name,obj.SimDist);
ReportTitle = sprintf('IRF Report:\\\\%s, %s',obj.Name,obj.SimDist);

%% -------------------------------------------------------------------

%% Prepare Draws
nDraws = obj.SimNDraws;
if strcmp(obj.SimDist,'PriorDraws')
    xd = obj.DrawPrior(nDraws);
elseif strcmp(obj.SimDist,'PostDraws')
    load(obj.FileName.MCMCDrawsRedux,'xd')
    xd = xd(:,randi(size(xd,2),1,nDraws));
else
    if ~isfield(obj.Param,obj.SimDist)
        fprintf(2,'Did not recognize distribution to use. Cannot proceed.\n');
        return
    end
    xd = obj.Param.(obj.SimDist);
end

%% Generate IRF
fprintf('Generating IRFs...\n');
% fnmats = @(x)feval(obj.FileName.Mats,x,...
%                'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
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
for jS=1:nShocks2Show
    Sj = obj.Shocks2Show{jS};
    for jP = 1:nPanels
        Pj = obj.FigPanels(jP);
        PlotData = nan(nDrawsUsed,nSteps,Pj.NVar);
        for jV=1:Pj.NVar
            Vj = Pj.Var{jV};
            [tf,idxV] = ismember(Vj,VarNames);
            if tf
               PlotData(:,:,jV) = Pj.Scale(jV)*obj.IRFShockSize(jS)*...
                   squeeze(IRF(idxV,:,jS,:))';
            end
        end
        PlotData = obj.VarScale*PlotData;
        Figj = Fig;
        Figj.TitleList = Pj.PrettyNames;
        if isfield(Pj,'FigShape');
            Figj.FigShape = Pj.FigShape;
        end
        OutFigj = vcFigure(PlotData,Figj);
        vcPrintPDF([obj.PlotDir.IRF,PlotFileName,...
             '_',Pj.Title,'_',Sj],Fig.KeepEPS,Fig.OpenPDF)
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
%         fprintf(fid,'\\includegraphics[width = \\textwidth,');
%         fprintf(fid,'viewport = 60 200 560 600, clip]{%s%s_%s_%s.pdf}\n',...
%                 obj.PlotDir.IRF,PlotFileName,Pj,Sj);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)

%% -------------------------------------------------------------------

if strcmp(Fig.Visible,'off')
    close all
end

%% Finish up
obj = obj.TrackTime(action,0);

end
%% -------------------------------------------------------------------
