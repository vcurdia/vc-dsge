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
op.Horizons = [1:32,inf];
op.Silent = 1;
op.VDPrctile = 50;
op.Fig.Visible = 'off';
op.Fig.Color = colorscheme('nColors',nShockVar,'LightFactors',[0,0.4,0.6]);
op.Fig.XTick = [1,4:4:32,33];
op.Fig.XTickLabels = {1,4:4:32,'inf'};
op.Fig.ShowPlotTitle=1;
op.Fig.LegPos = 'EO';
op.Fig.LegOrientation = 'vertical';
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

%% Prepare for IRF
if nargin<2 || isempty(xd)
    xd = obj.Param.Values;
end
nDraws = size(xd,2);

%% Generate VD
fnmats = @(x)obj.mats(x,...
               'StoreParam',0,'StoreStateEq',0,'StoreKF',0,'StoreAuxEq',0);
VDCheck = ones(1,nDraws);
nStateVar = obj.StateVar.N;
nObsVar = obj.ObsVar.N;
nAuxVar = obj.AuxVar.N;
nShockVar = obj.ShockVar.N;
isInfHorizon = ismember(inf,op.VDHorizons);
VDHorizons = sort(op.VDHorizons);
nHorizons = length(VDHorizons);
MaxHorizon = VDHorizons(end-isInfHorizon);
idxMat = eye(nShockVar);
VD = nan(nStateVar+nObsVar+nAuxVar,nHorizons,nShockVar,nDraws);
% ShockIdx = zeros(nShocks2Show,1);
% for j = 1:nShockVar
%     ShockIdx(j) = find(ismember(obj.ShockVar.Names,op.Shocks2Show(j)));
% end   
vd = cell(1,nDraws);
VDCheck = cell(1,nDraws);
parfor jd=1:nDraws
    matj = fnmats(xd(:,jd));
    checkj = all(matj.REE.eu==1);
    if checkj
        Vj = zeros(nStateVar+nObsVar+nAuxVar,nHorizons,nShockVar);
        for jS=1:nShockVar
            V = matj.REE.G2*idxMat(:,jS)*idxMat(jS,:)*matj.REE.G2';
            for jH=2:MaxHorizon
                V(:,:,jH) = matj.REE.G1*Vs(:,:,jH-1)*matj.REE.G1'+V(:,:,1);
            end
            V = V(:,:,VDHorizons(1:end-isInfHorizon));
            if isInfHorizon
                V(:,:,nHorizons) = real(...
                    lyapcsdSilent(matj.REE.G1,Vs(:,:,1),op.Silent));
            end
            if nObsVar>0
                for jH=1:nHorizons
                    VObs(:,:,jH) = matj.ObsEq.H*V(:,:,jH)*matj.ObsEq.H';
                end
            end
            if nAuxVar>0
                for jH=1:nHorizons
                    VAux(:,:,jH) = matj.AuxREE.G1*V(:,:,jH)*matj.AuxREE.G1';
                end
            end
            for jH=1:nHorizons
                v = diag(V(:,:,jH));
                if nObsVar>0, v = [v;diag(VObs(:,:,jH))]; end
                if nAuxVar>0, v = [v;diag(VAux(:,:,jH))]; end
                Vj(:,jS,jH) = v;
            end
        end
        vd{jd} = abs(Vj./repmat(sum(Vj,3),[1,1,nShockVar]));
        vcCheck{jd} = 1;
    else
        vdCheck{jd} = 0;
    end
end
for jd=1:nDraws, VD(:,:,:,jd) = vd{jd}; end clear vd
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
if nDrawsUsed>1
    fprintf('Percentile: %.1f\n',op.VDPrctile)
end
fprintf('\n')

HERE

for jPanel=1:nPanels
    fprintf('\nPanel: %s\n',Panels{jPanel})
    Vars2Showj = Vars2Show{jPanel};
    nVars2Showj = length(Vars2Showj);
    ShowVDj = ShowVD{jPanel};
    znamelength = [cellfun('length',Vars2Showj)];
    znamelengthmax = max(znamelength);
    snamelength = [cellfun('length',ShockVar)];
    snamelengthmax = max(max(snamelength),5);
    for jh=xTickIdx
        if nDrawsUsed==1
            VDj = ShowVDj(:,:,jh);
        else
            VDj = prctile(ShowVDj(:,:,jh,:),VDPrctile,4);
        end
        fprintf('\nHorizon: %.0f\n',VDHorizons(jh))
        fprintf(['%-',int2str(znamelengthmax),'s',...
                 repmat(['   %',int2str(snamelengthmax),'s'],1,nShockVar),...
                 '\n'],'',ShockVar{:})
        for jz=1:nVars2Showj
            fprintf(['%-',int2str(znamelengthmax),'s',...
                     repmat(['   %',int2str(snamelengthmax),'.3f'],1,...
                            nShockVar),'\n'],Vars2Showj{jz},VDj(jz,:))
        end
    end
end



%% Plot VDs
fprintf('Plotting VDs...\n');
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
