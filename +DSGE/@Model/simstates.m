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
if obj.AuxVar.N>0
    op.FigPanels = obj.setvarfigpanels('PanelList',{'StateVar','AuxVar'});
else
    op.FigPanels = obj.setvarfigpanels('PanelList',{'StateVar'});
end
op.TimeIdx = data.TimeIdx;
op.Tick = data.Tick;
op.TickLabel = data.TickLabel;
op.Fig.Visible = 'off';
op.PlotDir = 'Plots_States/';

op = updateoptions(op,varargin{:});


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
StatesCheck = zeros(1,nDraws);
parfor jd=1:nDraws
    mats = obj.mats(xd(:,jd));
    checkj = all(mats.REE.eu==1);
    if ~checkj
        StatesCheck(jd) = 0;
        continue
    end
    dj = dksmoother(mats,data,op.DrawStates);
    sj = [dj.StateVar;
          mats.KF.ObsVarBar+mats.ObsEq.H*dj.StateVar];
    if nAuxVar>0
        sj = [sj;
              mats.AuxREE.GBar ...
              + mats.AuxREE.G1*[dj.StateVar0,dj.StateVar(:,1:T-1)] ...
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
VarNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
nPanels = length(op.FigPanels);
TimeIdx = op.TimeIdx(ismember(op.TimeIdx,data.TimeIdx));
tid = ismember(data.TimeIdx,TimeIdx);
ntid = length(TimeIdx);
Fig.XTick = find(ismember(op.Tick,TimeIdx);
Fig.XTickLabel = op.TickLabel(ismember(TimeIdx,op.TickLabel));
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
    PlotData = nan(nDrawsUsed,TimeIdx,nVar);
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
        print('-dpdf',[op.PlotDir,PlotFileName,'_',Pj.Title])
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
%         fprintf(fid,'\\includegraphics[width=\\textwidth]{%s%s_%s_%s.pdf}\n',...
%                 op.PlotDir,PlotFileName,Pj,Sj);
        fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
                     '130 230 490 560]{%s%s_%s_%s.pdf}\n'],...
                op.PlotDir,PlotFileName,Pj,Sj);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end
fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
obj.TimeElapsed.stop(ttName)
