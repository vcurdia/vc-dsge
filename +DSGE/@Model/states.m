function states(obj,xd,varargin)

% states
% 
% Simulate states of DSGE model
% 
% See also:
% DSGE, SetupMyDSGE
%
% .............................................................................
% 
% Created: April 14, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Options
op.FNSuffix = '';
op.Data = obj.Data;
op.DrawStates = [];
op.Time2Show = [];
op.Tick.Labels = obj.Data.TickLabels;
op.Fig.Visible = 'off';
% op.Fig.YMinScale = 0.01;
op.Fig.Plot.LineWidth = 1.5;
op.PlotDir = 'plots-states/';
op.FigPanelsOptions = struct;
op.FigPanelsOptions.Shape = {3,1};

op = updateoptions(op,varargin{:});

if ~isempty(op.Data)
    data = op.Data;
else
    fprintf('Data is empty. Cannot simulate states.\n')
    return
end
if isempty(op.Time2Show), op.Time2Show = data.TimeIdx([1,end]); end

if ~isfield(op,'FigPanels')
    panelList = {'StateVar'};
    if obj.ObsVar.N>0
        panelList = [{'ObsVar'},panelList];
    end
    if obj.AuxVar.N>0
        panelList = [panelList,{'AuxVar'}];
    end
    op.FigPanels = obj.figpanels(op.FigPanelsOptions,'PanelList',panelList);
end


%% Preamble

fprintf('\nSimulating States\n')

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s-states%s',obj.Name,op.FNSuffix); 
ReportFileName = sprintf('report-%s-states%s',obj.Name,op.FNSuffix);
ReportTitle = sprintf('%s\\\\[30pt]States\\\\%s',obj.Name,...
                      strrep(op.FNSuffix,'-',''));

if nargin<2 || isempty(data)
    error('Cannot simulate states without data')
end

if nargin<3 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

if isempty(op.DrawStates),op.DrawStates = (nDraws>1); end

nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;


%% simulate states
States = nan(nStateVar+nObsVar+nAuxVar,data.T,nDraws);
StatesCheck = ones(1,nDraws);
parfor jd=1:nDraws
    mats = obj.solveree(xd(:,jd));
    checkj = mats.Status;
    if ~checkj
        StatesCheck(jd) = 0;
        continue
    end
    dj = dksmoother(mats,data.Values,op.DrawStates);
    dj.StateVar = mats.KF.StateVarBar + dj.StateVar;
    dj.StateVar0 = mats.KF.StateVarBar + dj.StateVar0;
    sj = [dj.StateVar; mats.ObsEq.HBar+mats.ObsEq.H*dj.StateVar];
    if nAuxVar>0
        sj = [sj; mats.AuxEq.PhiBar + mats.AuxEq.Phi*dj.StateVar];
    end
    States(:,:,jd) = sj;
end
States(:,:,~StatesCheck) = [];
StatesCheck(~StatesCheck) = [];
nDrawsUsed = length(StatesCheck);


%% Plot States
fprintf('Plotting States...\n');
Fig = op.Fig;
Fig.PlotBands = (nDraws>1);
vNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
nPanels = length(op.FigPanels);
tid = timeidx(op.Time2Show{:});
tid = tid(ismember(tid,data.TimeIdx));
idxT = ismember(data.TimeIdx,tid);
T = length(tid);
[Fig.XTick,Fig.XTickLabel] = setticklabel(tid,op.Tick);
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
    PlotData = nan(nDrawsUsed,T,nVar);
    for jV=1:nVar
        Vj = Pj.Names{jV};
        [tf,idxV] = ismember(Vj,vNames);
        if tf
            PlotData(:,:,jV) = Pj.Scale(jV)*squeeze(States(idxV,idxT,:))';
        end
    end
    h = vcfigure(PlotData,Figj);
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
