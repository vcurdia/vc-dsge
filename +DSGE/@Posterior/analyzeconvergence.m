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
op.Draws.BurnIn = 0.5;
op.Draws.Thinning = 1;
op.Draws.AuxParam = 1;
op.Table = DSGE.Options.Table;
op.NBin = 50;
op.Fig.Visible = 'off';
op.Fig.Color = colorscheme;
op.Fig.FontSize = 6;
op.PlotDir = 'Plots_Convergence/';

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Analyzing convergence of MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('AnalyzeConvergenceMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
ReportFileName = sprintf('%s_Report_MCMC_%.0f_Convergence',obj.Model.Name,...
                         obj.MCMCStage);
ReportTitle = sprintf('%s\\\\MCMC Stage %.0f\\\\Convergence Analysis',...
                      obj.Model.Name,obj.MCMCStage);

sample = obj.MCMCSample;

%% load the mcmc draws
draws = obj.loaddraws(op.Draws,'BurnIn',0,'CombineChains',0,'ExpandParam',0);


%% Plot draws
pdFN = sprintf('%s%s_Plots_MCMC_%.0f_Draws',...
             op.PlotDir,obj.Model.Name,obj.MCMCStage);
pdList = {};
pdListPretty = {};
pdNames = {};
pdPrettyNames = {};

%% LPDF
pdList{1} = 'LPDF';
pdListPretty{1} = 'Log-PDF';
pdNames{1} = pdList(1);
pdPrettyNames = pdListPretty(1);
pMax = max(draws.LPDF(:));
pMin = min(draws.LPDF(:));
pSpread = pMax-pMin;
pMax = pMax+.01*pSpread;
pMin = pMin-.01*pSpread;
vcfigure(draws.LPDF,op.Fig,'Shape',[sample.NChains,1],'YLim',[pMin,pMax]);
subplot(sample.NChains,1,1)
title(sprintf('Posterior %s in each chain',pdListPretty{1}))
printpdf(sprintf('%s_%s',pdFN,pdList{1}))


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
    fprintf(fid,'\\subsection{%s}\n',pdListPretty{jL});
    nFig = length(pdNames{jL});
    for jF=1:nFig
        if nFig>1
            fprintf(fid,'\\subsection{%s}\n',pdPrettyNames{jL}{jF});
        end
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{Fig_%s}\n',pdNames{jL}{jF});
        fprintf(fid,'\\includegraphics[width=\\textwidth]{%s_%s.pdf}\n',...
                pdFN,pdNames{jL}{jF});
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\clearpage \n');
    end
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)


%% Finish up
obj.TimeElapsed.stop(ttName)

