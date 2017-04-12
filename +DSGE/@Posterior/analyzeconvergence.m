function analyzeconvergence(obj,varargin)

% analyzeconvergence
% 
% Analyze convergence of MCMC sample
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 8, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.Draws.BurnIn = 0.25;
op.Draws.Thinning = 1;
op.Draws.AuxParam = 1;
op.Table = DSGE.Options.Table;
op.NBin = 50;
op.TraceStep = [];
op.Fig.Visible = 'off';
op.Fig.Color = colorscheme;
op.Fig.FontSize = 6;
op.PlotDirDraws = 'Plots_Draws/';
op.PlotDirTrace = 'Plots_Trace/';

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Analyzing convergence of MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('AnalyzeConvergenceMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

if ~isdir(op.PlotDirDraws),mkdir(op.PlotDirDraws),end
if ~isdir(op.PlotDirTrace),mkdir(op.PlotDirTrace),end
ReportFileName = sprintf('%s_Report_MCMC_%.0f_Convergence',obj.Model.Name,...
                         obj.MCMCStage);
ReportTitle = sprintf('%s\\\\MCMC Stage %.0f\\\\Convergence Analysis',...
                      obj.Model.Name,obj.MCMCStage);

sample = obj.MCMCSample;

%% load the mcmc draws
draws = obj.loaddraws(op.Draws,'BurnIn',0,'CombineChains',0,'ExpandParam',0);
nDrawsUsed = size(draws.LPDF,2);

if isempty(op.NBin), op.NBin = round(2*nDrawsUsed^(1/3)); end

p.LPDF.Title = 'Log-PDF';
p.LPDF.Names = {'LPDF'};
p.LPDF.PrettyNames = {p.LPDF.Title};
p.Param.Title = 'Parameters';
p.Param.Names = obj.Model.Param.Names(obj.EstimateIdx);
p.Param.PrettyNames = obj.Model.Param.PrettyNames(obj.EstimateIdx);
p.AuxParam.Title = 'Auxiliary Parameters';
p.AuxParam.Names = obj.Model.AuxParam.Names;
p.AuxParam.PrettyNames = obj.Model.AuxParam.PrettyNames;
pList = fieldnames(p);
nList = length(pList);
for jL=1:nList
    Lj = pList{jL};
    p.(Lj).N = length(p.(Lj).Names);
end


%% Plot draws
pdFN = sprintf('%s%s_Plots_MCMC_%.0f_Draws',...
             op.PlotDirDraws,obj.Model.Name,obj.MCMCStage);
pdList = pList;
for jL=1:length(pdList)
    Lj = pdList{jL};
    nc = 2 - strcmp(Lj,'LPDF');
    for jF=1:p.(Lj).N
        xd = draws.(Lj)(jF,:,:);
        pMax = max(xd(:));
        pMin = min(xd(:));
        pSpread = pMax-pMin;
        pMax = pMax+.01*pSpread;
        pMin = pMin-.01*pSpread;
        figure('Visible',op.Fig.Visible)
        clear h
        for jChain=1:sample.NChains
            h(jChain,1) = subplot(sample.NChains,nc,(jChain-1)*nc+1);
            plot(xd(1,:,jChain))
            ylim([pMin pMax])
            h(jChain,1).XTick = [op.Draws.BurnIn,1]*nDrawsUsed;
            h(jChain,1).XGrid = 'on';
            h(jChain,1).FontSize = op.Fig.FontSize;
            if nc>1
                h(jChain,2) = subplot(sample.NChains,2,jChain*2);
                histogram(xd(1,op.Draws.BurnIn*nDrawsUsed+1:end,jChain),...
                          op.NBin,'Normalization','probability',...
                          'FaceColor',op.Fig.Color(1,:),'FaceAlpha',1)
                xlim([pMin pMax])
                h(jChain,2).FontSize = op.Fig.FontSize;
            end
        end
        title(h(1,1),sprintf('%s in each chain',p.(Lj).PrettyNames{jF}))
        if nc>1
            title(h(1,2),sprintf('Hist excluding initial %.0f\\%% of obs',...
                                 100*op.Draws.BurnIn))
        end
        print('-dpdf',sprintf('%s_%s.pdf',pdFN,p.(Lj).Names{jF}))
    end
end
close all

%% eliminate BurnIn
draws.Param = draws.Param(:,(nDrawsUsed*op.Draws.BurnIn+1):end,:);
draws.LPDF = draws.LPDF(:,(nDrawsUsed*op.Draws.BurnIn+1):end,:);
draws.AuxParam = draws.AuxParam(:,(nDrawsUsed*op.Draws.BurnIn+1):end,:);
nDrawsUsed = size(draws.LPDF,2);

%% Plot trace
ptFN = sprintf('%s%s_Plots_MCMC_%.0f_Draws',...
             op.PlotDirTrace,obj.Model.Name,obj.MCMCStage);
ptList = {'Param','AuxParam'};
if isempty(op.TraceStep)
    op.TraceStep = floor(max(nDrawsUsed/100,1));
end
for jL=1:length(ptList)
    Lj = ptList{jL};
    for jF=1:p.(Lj).N
        xd = squeeze(draws.(Lj)(jF,:,:));
        SampleID = [op.TraceStep:op.TraceStep:nDrawsUsed]';
        nSample = length(SampleID);
        clear RollingMean RollingSD
        for js=1:nSample;
            RollingMean(js,:) = mean(xd(1:SampleID(js),:),1);
            RollingSD(js,:) = std(xd(1:SampleID(js),:),0,1);
        end
        pMean = mean(xd(:));
        pSD = std(xd(:));
        MeanBounds = [min(min(RollingMean(:)),pMean-2*pSD),...
                      max(max(RollingMean(:)),pMean+2*pSD)];
        MeanBounds = MeanBounds + ...
            (-1).^(1:-1:0)*0.01*(MeanBounds(2)-MeanBounds(1));
        SDBounds = [0 1.01*max(max(RollingSD(:)),2*pSD)];
        figure('Visible',op.Fig.Visible)
        clear h
        for jChain=1:sample.NChains
            h(jChain,1) = subplot(sample.NChains,2,(jChain-1)*nc+1);
            plot(SampleID,pMean*ones(size(SampleID)),'-',...
                 'Color',op.Fig.Color(1,:))
            hold on
            plot(SampleID,(pMean-pSD)*ones(size(SampleID)),':',...
                 'Color',op.Fig.Color(2,:))
            plot(SampleID,(pMean+pSD)*ones(size(SampleID)),':',...
                 'Color',op.Fig.Color(2,:))
            plot(SampleID,RollingMean(:,jChain),'-','Color',op.Fig.Color(1,:),...
                 'LineWidth',2)
            ylim(MeanBounds)
            xlim(SampleID([1,end]))
            h(jChain,1).FontSize = op.Fig.FontSize;
            h(jChain,2) = subplot(sample.NChains,2,jChain*2);
            plot(SampleID,pSD*ones(size(SampleID)),'-',...
                 'Color',op.Fig.Color(1,:))
            hold on
            plot(SampleID,RollingSD(:,jChain),'-','Color',op.Fig.Color(1,:),...
                 'LineWidth',2)
            ylim(SDBounds)
            xlim(SampleID([1,end]))
            h(jChain,2).FontSize = op.Fig.FontSize;
        end
        title(h(1,1),sprintf('Rolling Mean of %s',p.(Lj).PrettyNames{jF}))
        title(h(1,2),sprintf('Rolling SD of %s',p.(Lj).PrettyNames{jF}))
        print('-dpdf',sprintf('%s_%s.pdf',ptFN,p.(Lj).Names{jF}))
    end
end
close all



%% create report
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);

fprintf(fid,'\\begin{equation*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'number of chains: & %.0f\\\\\n',sample.NChains);
fprintf(fid,'size of each chain: & %.0f\\\\\n',sample.NDraws);
fprintf(fid,'burn in used: & %.0f (%.0f\\%%)\\\\\n',...
        op.Draws.BurnIn*sample.NDraws,op.Draws.BurnIn*100);
fprintf(fid,'thinning used: & %.0f\\\\\n',op.Draws.Thinning);
fprintf(fid,'number of draws used: & %.0f\\\\\\\\\n',draws.N);
fprintf(fid,'log-marginal likelihood: & %.4f\n',obj.LogMgLikelihood);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{equation*}\n');
fprintf(fid,'\\newpage\n');

fprintf(fid,'\\section{Plot Draws}\n');
for jL=1:length(pdList)
    Lj = pdList{jL};
    fprintf(fid,'\\subsection{%s}\n',p.(Lj).Title);
    for jF=1:p.(Lj).N
        if p.(Lj).N>1
            fprintf(fid,'\\subsubsection{%s}\n',p.(Lj).Names{jF});
        end
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{Fig_%s}\n',p.(Lj).Names{jF});
        fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
                     '130 230 490 540]{%s_%s.pdf}\n'],pdFN,p.(Lj).Names{jF});
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\clearpage \n');
    end
end

fprintf(fid,'\\section{Trace Plots}\n');
for jL=1:length(ptList)
    Lj = ptList{jL};
    fprintf(fid,'\\subsection{%s}\n',p.(Lj).Title);
    for jF=1:p.(Lj).N
        if p.(Lj).N>1
            fprintf(fid,'\\subsubsection{%s}\n',p.(Lj).Names{jF});
        end
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{Fig_%s}\n',p.(Lj).Names{jF});
        fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
                     '130 230 490 540]{%s_%s.pdf}\n'],ptFN,p.(Lj).Names{jF});
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\clearpage \n');
    end
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)

%% Finish up
obj.TimeElapsed.stop(ttName)

