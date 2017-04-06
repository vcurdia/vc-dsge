function analyzeparam(obj,varargin)

% analyzeparam
% 
% Analyze parameter posterior MCMC sample
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.BurnIn = 0.5;
op.Thinning = 1;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op.Table = DSGE.Options.Table;
op.NBin = 50;
op.Fig = DSGE.Options.Figure;
op.FigShape = [3,3];
op.LineColor = colorscheme;
op.PlotDir = 'Plots/PriorPost/';

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Analyzing MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('AnalyzeParamMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

ReportFileName = sprintf('%s_Report_MCMC_%.0f_Param',obj.Model.Name,...
                         obj.MCMCStage);
ReportTitle = sprintf('%s\\\\MCMC Stage %.0f\\\\Parameter Analysis',...
                      obj.Model.Name,obj.MCMCStage);

sample = obj.MCMCSample;
np = obj.Model.Param.N;
pNames = obj.Model.Param.Names;
nAux = obj.Model.AuxParam.N;
auxNames = obj.Model.AuxParam.Names;

%% load the mcmc draws
idxDraws = (op.BurnIn*sample.NDraws+1):op.Thinning:sample.NDraws;
for jChain=1:sample.NChains
    load(sample.FileNameDraws{jChain})
    xd(:,:,jChain) = xDraws(:,idxDraws);
    lpdfd(:,:,jChain) = lpdfDraws(:,idxDraws);
end
nDrawsUsed = size(xd,2)*sample.NChains;
xd = obj.expandparam(reshape(xd,obj.NEstimate,nDrawsUsed));
lpdfd = reshape(lpdfd,1,nDrawsUsed);
fprintf('Total number of draws per chain: %.0f\n', sample.NDraws)
fprintf('Burn in: %.0f%%\n', 100*op.BurnIn)
fprintf('Thinning used: %.0f\n', op.Thinning)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% Generate AuxParam
fprintf('Generating AuxParam draws...\n')
xdAux = zeros(nAux,nDrawsUsed);
fh = @(x)obj.Model.mats(x);
parfor jd=1:nDrawsUsed
    Matsj = fh(xd(:,jd));
    xdAux(:,jd) = Matsj.AuxParam;
end


%% Analyze sample
obj.MCMCSample.Param = sumstats(xd,op.Percentiles);
obj.MCMCSample.AuxParam = sumstats(xdAux,op.Percentiles);
Var = xd-repmat(obj.Mean,1,nDrawsUsed);
CorrMat = zeros(np);
for jr=1:np
    for jc=1:np
        CorrMat(jr,jc) = Var(jr,jc)/sqrt(Var(jr,jr)*Var(jc,jc));
    end
end
obj.Mean = obj.MCMCSample.Param.Mean;
obj.Median = obj.MCMCSample.Param.Median;
obj.SD = obj.MCMCSample.Param.SD;
obj.Prc05 = obj.MCMCSample.Param.Prc05;
obj.Prc95 = obj.MCMCSample.Param.Prc95;
obj.Var = Var*Var'/(nDrawsUsed-1);
obj.Corr = CorrMat;


%% check for new mode in mcmc sample
[lpdfNewMode,idxMax] = max(lpdfd);
fprintf('Checking MCMC draws for new mode\n')
fprintf('Previous mode lpdf: %.6f\n',obj.LPDFMode)
fprintf('Highest posterior density in MCMC draws: %.6f\n',lpdfNewMode)
if lpdfNewMode>obj.LPDFMode
    fprintf('Found MCMC draw with higher posterior density.\n')
    fprintf('Posterior mode updated.\n')
    namelength = [cellfun('length',pNames)];
    namelengthmax = max(namelength);
    fprintf(['%',int2str(namelengthmax),'s %-7s %-7s\n'],'','Old','New')
    for jp=1:np
        fprintf(['%',int2str(namelengthmax),'s %7.4f %7.4f\n'],...
                pNames{jp},obj.Mode(jp),xd(jp,idxMax))
    end
    obj.LPDFMode = lpdfNewMode;
    obj.Mode = xd(:,idxMax);
else
    fprintf('Did not find MCMC draw with higher posterior density.\n')
    fprintf('Previous posterior mode kept.\n')
end


%% Marginal likelihood
fprintf('Computing marginal likelihood...\n')
tau = 0.1:0.1:0.9;
ntau = length(tau);
pIdx = obj.EstimateIdx;
npd = obj.NEstimate;
% InvVar = inv(Post.Var);
xdd = xd(pIdx,:)-repmat(obj.Mean(pIdx),1,nDrawsUsed);
xddvar = xdd*xdd'/nDrawsUsed;
[xddvaru,xddvars,xddvarv] = svd(xddvar);
xddvarmd = min(size(xddvars));
bigev = find(diag(xddvars(1:xddvarmd,1:xddvarmd))>1e-6);
xddvardim = length(bigev);
% [rank(xdvar),np,xdvardim]
xddvarlndet = 0;
for j=1:npd
    if j>xddvardim
        xddvars(j,j) = 0;
    else
        xddvarlndet = xddvarlndet+log(xddvars(j,j));
        xddvars(j,j) = 1/xddvars(j,j);
    end
end
InvVar = xddvaru*xddvars*xddvaru';
lpdfMax = obj.LPDFMode;
lpdfMean = mean(lpdfd);
% [postMax,postMean]
% Constant terms
% lfConst = -log(tau)-np/2*log(2*pi)-1/2*log(det(Post.Var));
lfConst = -log(tau)-npd/2*log(2*pi)-1/2*xddvarlndet;
chi2Crit = chi2inv(tau,npd);
% Calculate the ratio of f(x)/post(x)
pw = zeros(ntau,nDrawsUsed);
nB = 10;
for jB=1:nB
%     fprintf('Computing %.0f of %.0f...\n',jB,nB)
    for jd=1:nDrawsUsed/nB
        idx = nDrawsUsed/nB*(jB-1)+jd;
        pratio = xdd(:,idx)'*InvVar*xdd(:,idx);
        pw(:,idx) = (pratio<=chi2Crit).*...
            exp(lfConst-1/2*pratio-lpdfd(idx)+lpdfMean);
    end
end
% harmonic mean
% LogMgLikelihood = postMean-log(mean(pw,2));
LogMgLikelihood = zeros(ntau,1);
for jtau=1:ntau
  pwj = pw(jtau,:);
  pwj = pwj(~isinf(pwj));
  pwj = pwj(~isnan(pwj));
  LogMgLikelihood(jtau) = lpdfMean-log(mean(pwj,2));
end
obj.LogMgLikelihood = LogMgLikelihood(tau==0.5);
% obj.LogMgLikelihood = mean(LogMgLikelihood);
for jt = 1:ntau
    obj.MCMCSample.LogMgLikelihoodtau(jt).Tau = tau(jt);
    obj.MCMCSample.LogMgLikelihoodtau(jt).Value = LogMgLikelihood(jt);
end


%% display results on screen
fprintf('\nResults from MCMC sample analysis:')
fprintf('\n==================================\n')
pNames = obj.Model.Param.Names;
pNameLength = [cellfun('length',pNames)];
pNameLengthMax = max(pNameLength);
DispList = {'','',pNames;
            'Prior',' Dist',obj.Prior.Dist;
%             '','   Mode',obj.Prior.Mode;
            '','   Mean',obj.Prior.Mean;
            '','     SD',obj.Prior.SD;
            '','     5%',obj.Prior.Prc05;
            '',' Median',obj.Prior.Median;
            '','    95%',obj.Prior.Prc95;
            'Posterior','   Mode',obj.Mode;
            '','   Mean',obj.Mean;
            '','     SD',obj.SD;
            '','     5%',obj.Prc05;
            '',' Median',obj.Median;
            '','    95%',obj.Prc95;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(pNameLengthMax),'s'],DispList{1,jr});
    str2show = sprintf('%s  %-5s',str2show,DispList{2,jr});
    for jc=3:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:obj.Model.Param.N
    str2show = sprintf(['%',int2str(pNameLengthMax),'s'],DispList{1,3}{jp});
    str2show = sprintf('%s  %5s',str2show,DispList{2,3}{jp});
    for jc=3:nc
        str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\n')

auxN = obj.Model.AuxParam.N;
auxNames = obj.Model.AuxParam.Names;
auxNameLength = [cellfun('length',auxNames)];
auxNameLengthMax = max(auxNameLength);
DispList = {'','',auxNames;
            'Prior','   Mean',obj.Prior.Sample.AuxParam.Mean;
            '','     SD',obj.Prior.Sample.AuxParam.SD;
            '','     5%',obj.Prior.Sample.AuxParam.Prc05;
            '',' Median',obj.Prior.Sample.AuxParam.Median;
            '','    95%',obj.Prior.Sample.AuxParam.Prc95;
            'Posterior','   Mean',obj.MCMCSample.AuxParam.Mean;
            '','     SD',obj.MCMCSample.AuxParam.SD;
            '','     5%',obj.MCMCSample.AuxParam.Prc05;
            '',' Median',obj.MCMCSample.AuxParam.Median;
            '','    95%',obj.MCMCSample.AuxParam.Prc95;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(auxNameLengthMax),'s'],DispList{1,jr});
    for jc=2:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:auxN
    str2show = sprintf(['%',int2str(auxNameLengthMax),'s'],DispList{1,3}{jp});
    for jc=2:nc
        str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\n')

% show mg Likelihood
fprintf('\nMarginal likelihood:')
fprintf('\n====================\n')
for jt=1:ntau
    fprintf('tau = %3.1f, log-marginal likelihood = %9.4f\n',...
            tau(jt),LogMgLikelihood(jt))
end
disp(' ')

% show correlation matrix
fprintf('\nCorrelation matrix:')
fprintf('\n===================\n')
namelength = [cellfun('length',pNames)];
namelengthmax = max(namelength);
str2show = sprintf(['%',int2str(namelengthmax),'s'],'');
for jc=1:np
    str2show = sprintf(['%s  %',int2str(namelengthmax),'s'],str2show,...
                       pNames{jc});
end
disp(str2show)
for jr=1:np
    str2show = sprintf(['%-',int2str(namelengthmax),'s'],pNames{jr});
    for jc=1:np
        str2show = sprintf(['%s  %',int2str(namelengthmax),'.4f'],str2show,...
                           obj.Corr(jr,jc));
    end
    disp(str2show)
end
disp(' ')


%% Make Prior Post Plots
fprintf('Making Prior-Posterior Plots...')
fn = sprintf('%s_Plots_MCMC_%.0f_PriorPost',obj.Model.Name,obj.MCMCStage);
nPlots = prod(op.FigShape);
nFig = ceil(np/nPlots);
for jF=1:nFig
    figure('Visible',op.Fig.Visible)
    for jf=1:nPlots 
        jp = (jF-1)*nPlots+jf;
        if jp>np, break, end
        subplot(op.FigShape(1),op.FigShape(2),jf)
        xPost = xd(jp,:);
%         [xFreq,xOut] = hist(xPost,nBin);
        hf = histogram(xPost,op.NBin);
%         xStep = xOut(2)-xOut(1);
%         xOutMin = min(xOut);
%         xOutMax = max(xOut);
        xStep = hf.BinWidth;
        xOut = hf.BinEdges;
        xCrit = [obj.Prior.Sample.Param.Prc01(jp),...
                 obj.Prior.Sample.Param.Prc99(jp)];
        xPlot = [sort(xOut(1)-xStep:-xStep:xCrit(1)),...
                 Out,...
                 xOut(end)+xStep:xStep:xCrit(2)];
        nPlot = zeros(size(xPlot));
        nIdx = ismember(xPlot,xOut);
        nPlot(nIdx) = hf.Values;
        for jx=1:length(xPlot)
%             xPriorPdf(jx) = eval(sprintf(Params(jp).priorpdfcmd,xPlot(jx)));
            xPriorPdf(jx) = obj.Prior.PDFCmd(xPlot(jx));
        end
%         nPlot = nPlot*(max(xPriorPdf)-min(xPriorPdf))/(max(nPlot)-min(nPlot));
        nPlot = nPlot/sum(nPlot*xStep); % normalize hist to have unit area
        HERE HERE
        bar(xPlot,nPlot)
        hold on
        plot(xPlot,xPriorPdf,'r','LineWidth',2)
        hold off
        title(Params(jp).prettyname)
        xBounds = xPlot([1,end]);
        xBounds = xBounds+(-1).^(1:-1:0)*0.01*(xBounds(2)-xBounds(1));
        xlim(xBounds)
        xTickLabels = xBounds(1):(xBounds(2)-xBounds(1))/8:xBounds(2);
        yBounds = max(max(nPlot),max(xPriorPdf));
        yBounds = [0,yBounds+0.01*yBounds];
        ylim(yBounds)
        set(gca,'YTick',[],'XTick',xTickLabels([2,5,8]),'FontSize',8)
        clear xPost xFreq xOut xStep xOutMin xOutMax xPlot nIdx xPriorPdf xBounds xTickLabels yBounds
    end
    vcPrintPDF(sprintf('%s%sFig%.0f',PlotDir.PriorPost,FileName.PlotsPriorPost,jF))
end


%% create report
fprintf('Making report: %s\n',ReportFileName);
fid = createtex(ReportFileName,ReportTitle);

fprintf(fid,'\\begin{equation*} \n');
fprintf(fid,'\\begin{tabular}{rl} \n');
fprintf(fid,'number of chains: & %.0f\\\\\n',sample.NChains);
fprintf(fid,'size of each chain: & %.0f\\\\\n',sample.NDraws);
fprintf(fid,'burn in used: & %.0f (%.0f\\%%)\\\\\n',...
        op.BurnIn*sample.NDraws,op.BurnIn*100);
fprintf(fid,'thinning used: & %.0f\\\\\n',op.Thinning);
fprintf(fid,'number of draws used: & %.0f\\\\\\\\\n',nDrawsUsed);
fprintf(fid,'log-marginal likelihood: & %.4f\n',obj.LogMgLikelihood);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{equation*}\n');

fprintf(fid,'\\newpage \n');
fprintf(fid,'\\section{Tables}\n');
fprintf(fid,'\\subsection{Parameters}\n');
np = obj.Model.Param.N;
str = [' & %.',int2str(op.Table.Precision),'f'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\subsection{Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    if op.Table.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+4+1+6));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{6}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mode & Mean & SD & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Prior.Dist{jr});
        fprintf(fid,str,obj.Prior.Prc05(jr));
        fprintf(fid,str,obj.Prior.Median(jr));
        fprintf(fid,str,obj.Prior.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Mode(jr));
        fprintf(fid,str,obj.Mean(jr));
        fprintf(fid,str,obj.SD(jr));
        fprintf(fid,str,obj.Prc05(jr));
        fprintf(fid,str,obj.Median(jr));
        fprintf(fid,str,obj.Prc95(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\subsection{Table: Auxiliary Parameters}\n');
tableBreaks = settablebreaks(auxN,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if jBreak>1
        fprintf(fid,'\\subsection{Auxiliary Parameters (Cont)}\n');
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('c',1,1+3+1+5));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{3}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{5}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mean & SD & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Model.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.MCMCSample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.MCMCSample.AuxParam.SD(jr));
        fprintf(fid,str,obj.MCMCSample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.MCMCSample.AuxParam.Median(jr));
        fprintf(fid,str,obj.MCMCSample.AuxParam.Prc95(jr));
        fprintf(fid,' \\\\\n');
        if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
            fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
        end        
    end
    fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\end{equation*}\n');
    fprintf(fid,'\\clearpage\n');
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)



%% Finish up
obj.TimeElapsed.stop(ttName)

