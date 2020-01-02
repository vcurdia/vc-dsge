function vd(obj,xd,varargin)

% vd
% 
% Generates variance decomposition for DSGE model
% 
% See also:
% DSGE, SetupMyDSGE
%
% .............................................................................
% 
% Created: April 12, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Default Options
op.FNSuffix = '';
op.NSteps = 25;
op.TickStep = 4;
op.FigPanels = obj.figpanels;
op.VDHorizons = [1:34,inf];
op.Verbose = 0;
op.VDPrctiles = [50]; % can create many at the same time
op.Table = DSGE.Options.Table;
op.Fig.Visible = 'off';
op.Fig.Color = [];
op.Fig.XTick = [1,4,8:8:32,35];
op.Fig.XTickLabel = {1,4:4:32,'   inf'};
op.Fig.ShowPlotTitle = 1;
op.Fig.LegPos = 'EO';
op.Fig.LegOrientation = 'vertical';
% op.Fig.FontSize = 8;
op.TightFig = 1;
op.TightFigOptions = struct;
op.PaperSize = [6.5, 6.5];
op.PaperPosition = [0, 0, 6.5, 6.5];
op.PlotDir = 'plots-vd/';

op = updateoptions(op,varargin{:});


%% Preamble

fprintf('\nMaking VD\n')

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s-vd%s',obj.Name,op.FNSuffix); 
ReportFileName = sprintf('report-%s-vd%s',obj.Name,op.FNSuffix);
ReportTitle = sprintf('%s\\\\[30pt]VD\\\\%s',obj.Name,...
                      strrep(op.FNSuffix,'-',''));

nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;
nShockVar = obj.ShockVar.N;
isInfHorizon = ismember(inf,op.VDHorizons);
VDHorizons = sort(op.VDHorizons);
nHorizons = length(VDHorizons);
MaxHorizon = VDHorizons(end-isInfHorizon);

vNames = [obj.StateVar.Names;obj.ObsVar.Names;obj.AuxVar.Names];
vNameLength = [cellfun('length',vNames)];
vNameLengthMax = max(vNameLength);
sNames = obj.ShockVar.Names;
sNameLength = [cellfun('length',sNames)];
sNameLengthMax = max(max(sNameLength),5);

tList = {'StateVar'};
if nObsVar>0, tList = {'ObsVar',tList{:}}; end
if nAuxVar>0, tList = {tList{:},'AuxVar'}; end

tid = 1:nHorizons;
if isempty(op.Fig.XTick), op.Fig.XTick = tid; end
if isempty(op.Fig.XTickLabel), op.Fig.XTickLabel = VDHorizons; end
if isempty(op.Fig.Color)
    op.Fig.Color = colorscheme('nColors',nShockVar,'LightFactors',[0,0.4, ...
                        0.6]);
end


%% Prepare draws
if nargin<2 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

%% Generate VD
% fnmats = @(x)obj.mats(x,...
%                'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
VDCheck = ones(1,nDraws);
idxMat = eye(nShockVar);
VD = nan(nStateVar+nObsVar+nAuxVar,nShockVar,nHorizons,nDraws);
VDCheck = zeros(1,nDraws);
verbose = op.Verbose;
parfor jd=1:nDraws
%     matj = fnmats(xd(:,jd));
    matj = obj.solveree(xd(:,jd));
    checkj = matj.Status;
    if checkj
        Vj = zeros(nStateVar+nObsVar+nAuxVar,nShockVar,nHorizons);
        for jS=1:nShockVar
            V = zeros(nStateVar,nStateVar,nHorizons);
            VAux = zeros(nAuxVar,nAuxVar,nHorizons);
            VObs = zeros(nObsVar,nObsVar,nHorizons);
            V(:,:,1) = matj.REE.G2*idxMat(:,jS)*idxMat(jS,:)*matj.REE.G2';
            for jH=2:MaxHorizon
                V(:,:,jH) = V(:,:,1) + matj.REE.G1*V(:,:,jH-1)*matj.REE.G1';
            end
            if isInfHorizon
                V(:,:,nHorizons) = real(...
                    lyapcsdvb(matj.REE.G1,V(:,:,1),verbose));
            end
            if nObsVar>0
                for jH=1:nHorizons
                    VObs(:,:,jH) = matj.ObsEq.H*V(:,:,jH)*matj.ObsEq.H';
                end
            end
            if nAuxVar>0
                for jH=1:nHorizons
                    VAux(:,:,jH) = matj.AuxEq.Phi*V(:,:,jH-1)*matj.AuxEq.Phi';
                end
            end
            for jH=1:nHorizons
                v = diag(V(:,:,jH));
                if nObsVar>0, v = [v;diag(VObs(:,:,jH))]; end
                if nAuxVar>0, v = [v;diag(VAux(:,:,jH))]; end
                Vj(:,jS,jH) = v;
            end
        end
        VD(:,:,:,jd) = abs(Vj./repmat(sum(Vj,2),[1,nShockVar,1]));
        VDCheck(jd) = 1;
    else
        VDCheck(jd) = 0;
    end
end
VD(:,:,:,~VDCheck) = [];
VDCheck(~VDCheck) = [];
nDrawsUsed = length(VDCheck);

%% Create tables
fprintf('\nVariance decomposition:')
fprintf('\n=======================\n')
if ~isempty(op.FNSuffix)
    fprintf('%s\n',op.FNSuffix)
end
fprintf('\n')
for jH=op.Fig.XTick
    for jL=1:length(tList)
        Lj = tList{jL};
        fprintf('Horizon: %.0f, %s',VDHorizons(jH),Lj)
        vj = obj.(Lj).Names;
        nj = obj.(Lj).N;
        vIdx = ismember(vNames,vj);
        for jPrc=1:max(length(op.VDPrctiles)*(nDrawsUsed>1),1)
            if nDrawsUsed==1
                VDj = VD(vIdx,:,jH);
                fprintf('\n')
            else
                Prcj = op.VDPrctiles(jPrc);
                fprintf(', Percentile %.1f',Prcj)
                VDj = prctile(VD(vIdx,:,jH,:),Prcj,4);
                fprintf('\n',Prcj)
            end
            fprintf(['%-',int2str(vNameLengthMax),'s',...
                     repmat(['   %',int2str(sNameLengthMax),'s'],1,nShockVar),...
                     '\n'],'',sNames{:})
            for jV=1:nj
                fprintf(['%-',int2str(vNameLengthMax),'s',...
                         repmat(['   %',int2str(sNameLengthMax),'.3f'],1,...
                                nShockVar),'\n'],vj{jV},VDj(jV,:))
            end
            fprintf('\n')
        end
    end
end

%% Plot VD
for jP=1:length(op.FigPanels)
    Pj = op.FigPanels(jP);
    Figj = op.Fig;
    Figj.TitleList = Pj.PrettyNames;
    Figj.Shape = Pj.Shape;
    hf = figure('Visible',Figj.Visible);
    clear ha
    nVar = length(Pj.Names);
    for jV=1:nVar
        Vj = Pj.Names{jV};
        ha(jV) = subplot(Figj.Shape{:},jV);
        vIdx = ismember(vNames,Vj);
        if nDrawsUsed==1
            PlotData = squeeze(VD(vIdx,:,:));
        else
            PlotData = squeeze(prctile(VD(vIdx,:,:,:),50,4));
        end
        bar(tid,permute(PlotData,[2,1]),'stacked','BarWidth',1,...
            'EdgeColor','none')
        axis tight
        colormap(Figj.Color)
        if Figj.ShowPlotTitle
            title(Pj.PrettyNames{jV});
        end
        ylim([0,1])
        ha(jV).XTick = Figj.XTick;
        ha(jV).XTickLabel = Figj.XTickLabel;
%         ha(jV).FontSize = Figj.FontSize;
    end
    if prod([Figj.Shape{:}])==1
        hl = legend(obj.ShockVar.PrettyNames,'Location',op.Fig.LegPos);
        hl.Orientation = op.Fig.LegOrientation;
        if strcmp(op.Fig.LegPos,'SO')
            hl.Position(2) = 0;
        end
    else
        hl = legend(obj.ShockVar.PrettyNames,'Location','S');
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


%% Make report with VD
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);
fprintf(fid,'\\newpage \n');

for jH=op.Fig.XTick
    fprintf(fid,'\n\\section{Horizon %.0f}\n',VDHorizons(jH));
    for jL=1:length(tList)
        Lj = tList{jL};
        vj = obj.(Lj).Names;
        nj = obj.(Lj).N;
        vIdx = ismember(vNames,vj);
        tableBreaks = settablebreaks(nj,op.Table.MaxRows);
        idxPar = 0;
        nBreaks = length(tableBreaks);
        for jBreak=1:nBreaks
            idxPar = (idxPar(end)+1):tableBreaks(jBreak);
            for jPrc=1:max(length(op.VDPrctiles)*(nDrawsUsed>1),1)
                fprintf(fid,'\\subsection{%s',Lj);
                if nDrawsUsed==1
                    VDj = VD(vIdx,:,jH);
                else
                    Prcj = op.VDPrctiles(jPrc);
                    fprintf(fid,', Prc %.0f',Prcj);
                    VDj = prctile(VD(vIdx,:,jH,:),Prcj,4);
                end
                if nBreaks>1, fprintf(fid,' (%.0f/%.0f)',jBreak,nBreaks); end
                fprintf(fid,'}\n');
                fprintf(fid,'\\begin{equation*}\n');
                if op.Table.MoveLeft
                    fprintf(fid,'\\hspace{-0.5in}\n');
                end
                fprintf(fid,'\\begin{tabular}{l%s}\n',repmat('r',1,nShockVar));
                fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
                fprintf(fid,' & %s', obj.ShockVar.PrettyNames{:});
                fprintf(fid,'\n\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
                for jr=idxPar
                    fprintf(fid,'%s',obj.(Lj).PrettyNames{jr});
                    fprintf(fid,' & $%.3f$',VDj(jr,:));
                    fprintf(fid,' \\\\\n');
                    if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
                        fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
                    end        
                end
                fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
                fprintf(fid,'\\end{tabular}\n');
                fprintf(fid,'\\end{equation*}\n');
                fprintf(fid,'\\clearpage \n');
            end
        end
    end
end

fprintf(fid,'\\section{Plots, Median}\n');
for jP=1:length(op.FigPanels)
    Pj = op.FigPanels(jP).Title;
    fprintf(fid,'\\subsection{%s}\n',strrep(Pj,'_',': '));
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{VD_%s}\n',Pj);
    fprintf(fid,'\\includegraphics[width=\\textwidth]{%s%s-%s.pdf}\n',...
            op.PlotDir,PlotFileName,Pj);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
end


fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
