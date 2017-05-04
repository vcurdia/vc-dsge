function irf(obj,xd,varargin)

% irf
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


%% Default Options
op.FNSuffix = '';
op.NSteps = 25;
op.TickStep = 4;
op.FigPanels = obj.figpanels;
op.Shocks2Show = obj.ShockVar.Names;
op.ShockSize = [];
op.Fig.Visible = 'off';
op.Fig.Plot.LineWidth = 1.5;
op.PlotDir = 'Plots_IRF/';

%% Update options
op = updateoptions(op,varargin{:});

%% Check options
nShocks2Show = length(op.Shocks2Show);
if isempty(op.ShockSize), op.ShockSize = ones(1,nShocks2Show); end


%% Preamble

fprintf('\n*** Making IRF\n')
ttName = ['IRF',op.FNSuffix];
obj.TimeTracker.start(ttName)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s_IRF%s',obj.Name,op.FNSuffix); 
ReportFileName = sprintf('%s_Report_IRF%s',obj.Name,op.FNSuffix);
ReportTitle = sprintf('%s\\\\IRF\\\\%s',obj.Name,...
                      strrep(op.FNSuffix,'_',''));

%% Prepare for IRF
if nargin<2 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

%% Generate IRF
fnmats = @(x)obj.mats(x,...
               'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
nSteps = op.NSteps;
nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;
ShockIdx = zeros(nShocks2Show,1);
for j = 1:nShocks2Show
    ShockIdx(j) = find(ismember(obj.ShockVar.Names,op.Shocks2Show(j)));
end   
IRF = nan(nStateVar+nObsVar+nAuxVar,nSteps,nShocks2Show,nDraws);
IRFCheck = ones(1,nDraws);
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
Fig = op.Fig;
Fig.PlotBands = (nDraws>1);
Fig.XTick = 1:op.TickStep:nSteps;
Fig.XTickLabel = 0:op.TickStep:(nSteps-1);
VarNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
nPanels = length(op.FigPanels);
for jP = 1:nPanels
    Pj = op.FigPanels(jP);
    Figj = Fig;
    if isfield(Pj,'PrettyNames')
        Figj.TitleList = Pj.PrettyNames;
    else
        Figj.TitleList = Pj.Names;
    end
    if isfield(Pj,'FigShape');
        Figj.Shape = Pj.FigShape;
    end
    nVar = length(Pj.Names);
    PlotData = nan(nDrawsUsed,nSteps,nVar,nShocks2Show);
    for jV=1:nVar
        Vj = Pj.Names{jV};
        [tf,idxV] = ismember(Vj,VarNames);
        if tf
            for jS=1:nShocks2Show
                PlotData(:,:,jV,jS) = Pj.Scale(jV)*...
                    op.ShockSize(jS)*squeeze(IRF(idxV,:,jS,:))';
            end
        end
    end
    for jS=1:nShocks2Show
        if jS==1
            h = vcfigure(PlotData(:,:,:,jS),Figj);
        else
            h = vcfigureupdate(h,PlotData(:,:,:,jS));
        end
%         printpdf([op.PlotDir,PlotFileName,...
%                     '_',Pj.Title,'_',op.Shocks2Show{jS}])
        print('-dpdf',[op.PlotDir,PlotFileName,...
                    '_',Pj.Title,'_',op.Shocks2Show{jS}])
    end
end

%% Make report with IRF
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
for jS=1:nShocks2Show
    Sj = op.Shocks2Show{jS};
    fprintf(fid,'\\section{Shock: %s}\n',strrep(Sj,'_',''));
    for jP = 1:nPanels
        Pj = op.FigPanels(jP).Title;
        fprintf(fid,'\\subsection{%s}\n',strrep(Pj,'_',': '));
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{IRF_%s_%s}\n',Pj,Sj);
        fprintf(fid,'\\includegraphics[width=\\textwidth]{%s%s_%s_%s.pdf}\n',...
                op.PlotDir,PlotFileName,Pj,Sj);
%         fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
%                      '130 230 490 560]{%s%s_%s_%s.pdf}\n'],...
%                 op.PlotDir,PlotFileName,Pj,Sj);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
obj.TimeTracker.stop(ttName)
