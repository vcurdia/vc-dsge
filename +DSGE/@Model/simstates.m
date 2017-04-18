function simstates(obj,data,xd,varargin)

% simstates
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


%% Default Options
op.FileNameSuffix = '';
op.DrawStates = [];
op.Time2Show = data.TimeIdx([1,end]);
op.Tick.Labels = [];
op.Fig.Visible = 'off';
op.Fig.YMinScale = 0.01;
op.PlotDir = 'Plots_States/';
op.FigPanelsOptions = struct;

op = updateoptions(op,varargin{:});

if ~isfield(op,'FigPanels')
    if obj.AuxVar.N>0
        op.FigPanels = obj.setvarfigpanels(op.FigPanelsOptions,...
                                           'PanelList',{'StateVar','AuxVar'});
    else
        op.FigPanels = obj.setvarfigpanels(op.FigPanelsOptions,...
                                           'PanelList',{'StateVar'});
    end
end


%% Preamble

fprintf('\n*** Simulating States\n')
ttName = ['States',op.FileNameSuffix];
obj.TimeElapsed.start(ttName)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s_States%s',obj.Name,op.FileNameSuffix); 
ReportFileName = sprintf('%s_Report_States%s',obj.Name,op.FileNameSuffix);
ReportTitle = sprintf('%s\\\\States Report %s',obj.Name,...
                      strrep(op.FileNameSuffix,'_',''));

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
    mats = obj.mats(xd(:,jd));
    checkj = all(mats.REE.eu==1);
    if ~checkj
        StatesCheck(jd) = 0;
        continue
    end
    dj = dksmoother(mats,data.Values,op.DrawStates);
    sj = [dj.StateVar;
          mats.KF.ObsVarBar+mats.ObsEq.H*dj.StateVar];
    if nAuxVar>0
        sj = [sj;
              mats.AuxREE.GBar ...
              + mats.AuxREE.G1*[dj.StateVar0,dj.StateVar(:,1:data.T-1)] ...
              + mats.AuxREE.G2*dj.ShockVar;
             ];
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
    if isfield(Pj,'FigShape');
        Figj.Shape = Pj.FigShape;
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
    print('-dpdf',[op.PlotDir,PlotFileName,'_',Pj.Title])
end

%% Make report with IRF
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');
for jP = 1:nPanels
    Pj = op.FigPanels(jP).Title;
    fprintf(fid,'\\section{%s}\n',strrep(Pj,'_',': '));
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{States_%s}\n',Pj);
    fprintf(fid,['\\includegraphics[width=\\textwidth]{%s%s_%s.pdf}\n'],...
            op.PlotDir,PlotFileName,Pj);
%     fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
%                  '120 230 490 560]{%s%s_%s.pdf}\n'],...
%             op.PlotDir,PlotFileName,Pj);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
obj.TimeElapsed.stop(ttName)
