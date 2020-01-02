function out = irf(obj,xd,varargin)

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
op.FigPanels = [];
op.Shocks2Show = obj.ShockVar.Names;
op.ShockSize = [];
op.Fig.Visible = 'off';
op.Fig.Plot.LineWidth = 1.5;
op.PlotDir = 'plots-irf/';
op.MakePlots = 1;

%% Update options
op = updateoptions(op,varargin{:});

%% Check options
if isempty(op.FigPanels), op.FigPanels = obj.figpanels; end
nShocks2Show = length(op.Shocks2Show);
if isempty(op.ShockSize), op.ShockSize = ones(1,nShocks2Show); end


%% Preamble

fprintf('\nMaking IRF\n')

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s-irf%s',obj.Name,op.FNSuffix); 
ReportFileName = ['report-',obj.Name,'-irf',op.FNSuffix];
ReportTitle = sprintf('%s\\\\[30pt]IRF\\\\%s',obj.Name,...
                      strrep(op.FNSuffix,'-',''));

out = struct;
out.Var = merge(obj.StateVar,obj.ObsVar,obj.AuxVar);
out.Shock = subset(obj.ShockVar,op.Shocks2Show);
out.NSteps = op.NSteps;

%% Prepare for IRF
if nargin<2 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

%% Generate IRF
% fnmats = @(x)obj.solveree(x);
nSteps = op.NSteps;
nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;
% [tfShock,idxShock] = ismember(op.Shocks2Show,obj.ShockVar.Names);
ShockIdx = zeros(nShocks2Show,1);
for j = 1:nShocks2Show
    ShockIdx(j) = find(ismember(obj.ShockVar.Names,op.Shocks2Show(j)));
end   
IRF = nan(nStateVar+nObsVar+nAuxVar,nSteps,nShocks2Show,nDraws);
IRFCheck = ones(1,nDraws);
parfor jd=1:nDraws
%     matj = fnmats(xd(:,jd));
    matj = obj.solveree(xd(:,jd));
    checkj = matj.Status;
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
        for t=1:nSteps
            irfAux(:,:,t) = matj.AuxEq.Phi*irf(:,:,t);
        end
        IRFj = [IRFj;irfAux];
    end
    IRF(:,:,:,jd) = permute(IRFj,[1,3,2]);
end
IRF(:,:,:,~IRFCheck) = [];
IRFCheck(~IRFCheck) = [];
nDrawsUsed = length(IRFCheck);
if nDrawsUsed==0
    error('No valid simulations were generated. Cannot proceed.')
end
out.IRF = IRF;
out.IRFCheck = IRFCheck;
out.NDraws = nDrawsUsed;

%% Plot IRFs
if op.MakePlots
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
        if isfield(Pj,'Shape');
            Figj.Shape = Pj.Shape;
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
            print('-dpdf',[op.PlotDir,PlotFileName,...
                           '-',Pj.Title,'-',op.Shocks2Show{jS}])
        end
    end
end

%% Make report with IRF
if op.MakePlots
    fprintf('Making report: %s\n',ReportFileName);
    fid = createtex(ReportFileName,ReportTitle);
    fprintf(fid,'\\newpage \n');
    for jS=1:nShocks2Show
        Sj = op.Shocks2Show{jS};
        fprintf(fid,'\\section{Shock: %s}\n',strrep(Sj,'_',''));
        for jP = 1:nPanels
            Pj = op.FigPanels(jP).Title;
            fprintf(fid,'\\subsection{%s}\n',Pj);
            fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
            fprintf(fid,'\\label{irf-%s-%s}\n',Pj,Sj);
            fprintf(fid,'\\includegraphics[width=\\textwidth]{%s%s-%s-%s.pdf}\n',...
                    op.PlotDir,PlotFileName,Pj,Sj);
            fprintf(fid,'\\end{figure}\n');
            fprintf(fid,'\\newpage \n');
        end
    end
    fprintf(fid,'\\end{document}\n');
    fclose(fid);
    pdflatex(ReportFileName)
end

%% Finish up
close all
