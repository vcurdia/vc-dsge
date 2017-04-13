function makevd(obj,xd,varargin)

% makevd
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
op.FileNameSuffix = '';
op.NSteps = 25;
op.TickStep = 4;
op.FigPanels = obj.setvarfigpanels;
op.VDHorizons = [1:34,inf];
op.Silent = 1;
op.VDPrctiles = [50,5,95];
op.Table = DSGE.Options.Table;
op.Fig.Visible = 'off';
op.Fig.Color = [];
op.Fig.XTick = [1,4:4:32,35];
op.Fig.XTickLabel = {1,4:4:32,'   inf'};
op.Fig.ShowPlotTitle = 1;
op.Fig.LegPos = 'EO';
op.Fig.LegOrientation = 'vertical';
op.Fig.FontSize = 8;
op.PlotDir = 'Plots_VD/';

op = updateoptions(op,varargin{:});


%% Preamble

fprintf('\n*** Making VD\n')
ttName = ['VD',op.FileNameSuffix];
obj.TimeElapsed.start(ttName)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
PlotFileName = sprintf('%s_VD%s',obj.Name,op.FileNameSuffix); 
ReportFileName = sprintf('%s_Report_VD%s',obj.Name,op.FileNameSuffix);
ReportTitle = sprintf('%s\\\\VD Report %s',obj.Name,...
                      strrep(op.FileNameSuffix,'_',''));

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
fnmats = @(x)obj.mats(x,...
               'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
VDCheck = ones(1,nDraws);
idxMat = eye(nShockVar);
VD = nan(nStateVar+nObsVar+nAuxVar,nShockVar,nHorizons,nDraws);
vd = cell(1,nDraws);
VDCheck = cell(1,nDraws);
isSilent = op.Silent;
parfor jd=1:nDraws
    matj = fnmats(xd(:,jd));
    checkj = all(matj.REE.eu==1);
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
            if nAuxVar>0
                VAux(:,:,1) = ...
                    matj.AuxREE.G2*idxMat(:,jS)*idxMat(jS,:)*matj.AuxREE.G2';
                for jH=2:MaxHorizon
                    VAux(:,:,jH) = VAux(:,:,1) + ...
                        matj.AuxREE.G1*V(:,:,jH-1)*matj.AuxREE.G1';
                end
            end
%             V = V(:,:,VDHorizons(1:end-isInfHorizon));
%             VAux = VAux(:,:,VDHorizons(1:end-isInfHorizon));
            if isInfHorizon
                V(:,:,nHorizons) = real(...
                    lyapcsdsilent(matj.REE.G1,V(:,:,1),isSilent));
                VAux(:,:,nHorizons) = VAux(:,:,1) + ...
                        matj.AuxREE.G1*V(:,:,nHorizons)*matj.AuxREE.G1';
            end
            if nObsVar>0
                for jH=1:nHorizons
                    VObs(:,:,jH) = matj.ObsEq.H*V(:,:,jH)*matj.ObsEq.H';
                end
            end
            for jH=1:nHorizons
                v = diag(V(:,:,jH));
                if nObsVar>0, v = [v;diag(VObs(:,:,jH))]; end
                if nAuxVar>0, v = [v;diag(VAux(:,:,jH))]; end
                Vj(:,jS,jH) = v;
            end
        end
        vd{jd} = abs(Vj./repmat(sum(Vj,2),[1,nShockVar,1]));
        VDCheck{jd} = 1;
    else
        VDCheck{jd} = 0;
    end
end
for jd=1:nDraws, VD(:,:,:,jd) = vd{jd}; end, clear vd
VDCheck = [VDCheck{:}];
VD(:,:,:,~VDCheck) = [];
VDCheck(~VDCheck) = [];
nDrawsUsed = length(VDCheck);

%% Create tables
fprintf('\nVariance decomposition:')
fprintf('\n=======================\n')
if ~isempty(op.FileNameSuffix)
    fprintf('%s\n',op.FileNameSuffix)
end
for jH=op.Fig.XTick
    for jL=1:length(tList)
        Lj = tList{jL};
        fprintf('\nHorizon: %.0f, %s\n',VDHorizons(jH),Lj)
        vj = obj.(Lj).Names;
        nj = obj.(Lj).N;
        vIdx = ismember(vNames,vj);
        for jPrc=1:max(length(op.VDPrctiles)*(nDrawsUsed>1),1)
            if nDrawsUsed==1
                VDj = VD(vIdx,:,jH);
            else
                Prcj = op.VDPrctiles(jPrc);
                fprintf('Percentile: %.1f\n',Prcj)
                VDj = prctile(VD(vIdx,:,jH,:),Prcj,4);
            end
            fprintf(['%-',int2str(vNameLengthMax),'s',...
                     repmat(['   %',int2str(sNameLengthMax),'s'],1,nShockVar),...
                     '\n'],'',sNames{:})
            for jV=1:nj
                fprintf(['%-',int2str(vNameLengthMax),'s',...
                         repmat(['   %',int2str(sNameLengthMax),'.3f'],1,...
                                nShockVar),'\n'],vj{jV},VDj(jV,:))
            end
        end
    end
end

%% Plot VD
for jP=1:length(op.FigPanels)
    Pj = op.FigPanels(jP);
    Figj = op.Fig;
    Figj.TitleList = Pj.PrettyNames;
    Figj.Shape = Pj.FigShape;
    figure('Visible',Figj.Visible)
    for jV=1:Pj.N
        Vj = Pj.Names{jV};
        hf(jV) = subplot(Figj.Shape{:},jV);
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
        ax = gca;
        ax.XTick = Figj.XTick;
        ax.XTickLabel = Figj.XTickLabel;
        ax.FontSize = Figj.FontSize;
    end
    hl = legend(obj.ShockVar.PrettyNames,'Location',op.Fig.LegPos);
    if prod([Figj.Shape{:}])==1
        hl.Orientation = op.Fig.LegOrientation;
        if strcmp(op.Fig.LegPos,'SO')
            hl.Position(2) = 0;
        end
    else
        hl.Orientation = 'horizontal';
        legPos = hl.Position;
        xL = hf((Figj.Shape{1}-1)*Figj.Shape{2}+1).Position;
        xR = hf((Figj.Shape{1}-1)*Figj.Shape{2}).Position;
        legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
        legPos(2) = 0;
        hl.Position = legPos;
    end
    print('-dpdf',[op.PlotDir,PlotFileName,'_',Pj.Title])
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
            if nBreaks==1
                fprintf(fid,'\\subsection{%s}\n',Lj);
            else
                fprintf(fid,'\\subsection{%s (%.0f/%.0f)}\n',Lj,...
                        jBreak,nBreaks);
            end
            for jPrc=1:max(length(op.VDPrctiles)*(nDrawsUsed>1),1)
                if nDrawsUsed==1
                    VDj = VD(vIdx,:,jH);
                else
                    Prcj = op.VDPrctiles(jPrc);
                    fprintf(fid,'\\subsubsection{Percentile: %.0f}\n',Prcj);
                    VDj = prctile(VD(vIdx,:,jH,:),Prcj,4);
                end
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

fprintf(fid,'\\section{VD Plots: Median}\n');
for jP=1:length(op.FigPanels)
    Pj = op.FigPanels(jP).Title;
    fprintf(fid,'\\subsection{%s}\n',strrep(Pj,'_',': '));
    fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
    fprintf(fid,'\\label{VD_%s}\n',Pj);
    fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
                 '130 230 490 560]{%s%s_%s.pdf}\n'],...
            op.PlotDir,PlotFileName,Pj);
    fprintf(fid,'\\end{figure}\n');
    fprintf(fid,'\\newpage \n');
end


fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
close all
obj.TimeElapsed.stop(ttName)
