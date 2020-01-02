function sd(obj,xd,varargin)

% sd
% 
% shock decomposition of DSGE model states
% 
% See also:
% DSGE, SetupMyDSGE
%
% .............................................................................
% 
% Created: April 28, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Options
op.FNSuffix = '';
op.Data = obj.Data;
op.DrawStates = [];
op.ShowOther = 1;
op.Time2Show = [];
op.Tick.Labels = obj.Data.TickLabels;
op.Fig.Visible = 'off';
op.Fig.Color = [];
% op.Fig.YMinScale = 0.01;
op.Fig.ShowPlotTitle = 1;
op.Fig.LegPos = 'EO';
op.Fig.LegOrientation = 'vertical';
% op.Fig.FontSize = 8;
op.TightFig = 1;
op.TightFigOptions = struct;
op.PaperSize = [6.5, 6.5];
op.PaperPosition = [0, 0, 6.5, 6.5];
op.PlotDir = 'plots-sd/';
op.FigPanelsOptions = struct;
op.FigPanelsOptions.Shape = {3,1};

op = updateoptions(op,varargin{:});

if ~isempty(op.Data)
    data = op.Data;
else
    fprintf('Data is empty. Cannot generate shock decomposition.\n')
    return
end
if isempty(op.Time2Show), op.Time2Show = data.TimeIdx([1,end]); end

if ~isfield(op,'ShockGroups')
    op.ShockGroups = cell(obj.ShockVar.N,2);
    for j=1:obj.ShockVar.N
        op.ShockGroups(j,:) = {obj.ShockVar.PrettyNames{j},...
                            obj.ShockVar.Names(j)};
    end
end
nGroups = length(op.ShockGroups);

if isempty(op.Fig.Color)
    op.Fig.Color = colorscheme('nColors',nGroups+op.ShowOther,...
                               'LightFactors',[0,0.4,0.6]);
end

if ~isfield(op,'FigPanels')
    op.FigPanels = obj.figpanels(op.FigPanelsOptions);
end


%% Preamble

fprintf('\nShock Decomposition\n')

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s-sd%s',obj.Name,op.FNSuffix); 
ReportFileName = sprintf('report-%s-sd%s',obj.Name,op.FNSuffix);
ReportTitle = sprintf('%s\\\\[30pt]Shock Decomposition\\\\%s',obj.Name,...
                      strrep(op.FNSuffix,'-',''));

if nargin<2 || isempty(data)
    error('Cannot perform shock decomposition without data')
end

if nargin<3 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

if isempty(op.DrawStates),op.DrawStates = (nDraws>1); end

nStateVar = obj.StateVar.N;
nShockVar = obj.ShockVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;

shockIdx = false(nGroups,nShockVar);
for jG=1:nGroups
    shockIdx(jG,:) = ismember(obj.ShockVar.Names,op.ShockGroups{jG,2});
end

%% shock decomposition
SD = nan(nStateVar+nObsVar+nAuxVar,data.T,nGroups+op.ShowOther,nDraws);
simCheck = ones(1,nDraws);
parfor jd=1:nDraws
    mats = obj.solveree(xd(:,jd));
    checkj = mats.Status;
    if ~checkj
        simCheck(jd) = 0;
        continue
    end
    dj = dksmoother(mats,data.Values,op.DrawStates);
    sdj = nan(nStateVar+nObsVar+nAuxVar,data.T,nGroups);
    addstates = eye(nStateVar);
    if nObsVar>0, addstates = [addstates;mats.ObsEq.H]; end
    if nAuxVar>0, addstates = [addstates;mats.AuxEq.Phi]; end
    for jG=1:nGroups
        ej = dj.ShockVar(shockIdx(jG,:),:);
        G2j = mats.REE.G2(:,shockIdx(jG,:));
        sj = zeros(nStateVar,data.T);
        sj(:,1) = G2j*ej(:,1);
        for t=2:data.T
            sj(:,t) = mats.REE.G1*sj(:,t-1) + G2j*ej(:,t);
        end
        sdj(:,:,jG) = addstates*sj;
    end
    if op.ShowOther
        sdj(:,:,nGroups+1) = addstates*dj.StateVar-sum(sdj,3);
    end
    SD(:,:,:,jd) = sdj;
end
SD(:,:,:,~simCheck) = [];
simCheck(~simCheck) = [];
nDrawsUsed = length(simCheck);



%% Plot States
fprintf('Plotting shock decomposition...\n');
Fig = op.Fig;
Fig.PlotBands = (nDraws>1);
vNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
nPanels = length(op.FigPanels);
tid = timeidx(op.Time2Show{:});
tid = tid(ismember(tid,data.TimeIdx));
idxT = ismember(data.TimeIdx,tid);
T = length(tid);
[Fig.XTick,Fig.XTickLabel] = setticklabel(tid,op.Tick);
GroupLabels = {op.ShockGroups{:,1}};
if op.ShowOther, GroupLabels{end+1} = 'Other'; end
for jP = 1:nPanels
    Pj = op.FigPanels(jP);
    Figj = Fig;
    Figj.TitleList = Pj.PrettyNames;
    Figj.Shape = Pj.Shape;
    hf = figure('Visible',Figj.Visible);
    nVar = length(Pj.Names);
    clear ha
    for jV=1:nVar
        Vj = Pj.Names{jV};
        ha(jV) = subplot(Figj.Shape{:},jV);
        vIdx = ismember(vNames,Vj);
        if nDrawsUsed==1
            PlotData = Pj.Scale(jV)*squeeze(SD(vIdx,idxT,:,:));
        else
            PlotData = Pj.Scale(jV)*squeeze(median(SD(vIdx,idxT,:,:),4));
        end
        bar(1:T,max(0,PlotData),'stacked','EdgeColor','none')
        hold on
        bar(1:T,min(0,PlotData),'stacked','EdgeColor','none')
        hold off
        axis tight
        colormap(Figj.Color)
        if Figj.ShowPlotTitle
            title(Pj.PrettyNames{jV});
        end
        ha(jV).XTick = Figj.XTick;
        ha(jV).XTickLabel = Figj.XTickLabel;
%         ha(jV).FontSize = Figj.FontSize;
    end
    if prod([Figj.Shape{:}])==1
        hl = legend(GroupLabels,'Location',op.Fig.LegPos);
        hl.Orientation = op.Fig.LegOrientation;
        if strcmp(op.Fig.LegPos,'SO')
            hl.Position(2) = 0;
        end
    else
        hl = legend(GroupLabels,'Location','S');
        hl.Orientation = 'horizontal';
        legPos = hl.Position;
        legPos(1) = 0.5-legPos(3)/2;
        legPos(2) = 0;
        hl.Position = legPos;
    end
    hf.PaperSize = op.PaperSize;
    hf.PaperPosition = op.PaperPosition;
    if op.TightFig
        tightfig(hf,Figj.Shape,ha,op.TightFigOptions)
    end
    print('-dpdf',[op.PlotDir,PlotFileName,'-',Pj.Title])
end

%% Make report 
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
for jP = 1:nPanels
    Pj = op.FigPanels(jP).Title;
    fprintf(fid,'\\section{%s}\n',strrep(Pj,'_',': '));
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{States_%s}\n',Pj);
    fprintf(fid,['\\includegraphics[width=\\textwidth]{%s%s-%s.pdf}\n'],...
            op.PlotDir,PlotFileName,Pj);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
