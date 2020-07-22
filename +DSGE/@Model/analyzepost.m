function analyzepost(obj,varargin)

% analyzepost
% 
% Analyze parameter posterior MCMC sample
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017-2018 Vasco Curdia


%% Options
op.Draws.BurnIn = 0.25;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op.Table = DSGE.Options.Table;
op.NDrawsPrior = 20000;
op.NBin = 50;
op.Fig.Visible = 'off';
op.Fig.Shape = {3,3};
op.Fig.Color = colorscheme;
op.Fig.FontSize = 8;
op.TightFig = 1;
op.TightFigOptions = struct;
op.PaperSize = [6.5, 6.5];
op.PaperPosition = [0, 0, 6.5, 6.5];
op.PlotDir = 'plots-priorpost/';
op.showcorr = 0;

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\nAnalyzing MCMC Sample %.0f\n',obj.Post.MCMCStage)

if ~isdir(op.PlotDir),mkdir(op.PlotDir),end
ReportFileName = sprintf('report-%s-param-post-mcmc-%.0f',obj.Name,...
                         obj.Post.MCMCStage);
ReportTitle = sprintf('%s\\\\[30pt]Parameter Analysis\\\\Posterior MCMC Stage %.0f',...
                      obj.Name,obj.Post.MCMCStage);

sample = obj.Post.MCMCSample;
np = obj.Param.N;
pNames = obj.Param.Names;
nAux = obj.AuxParam.N;
auxNames = obj.AuxParam.Names;

%% load the mcmc draws
op.Draws.AuxParam = 1;
draws = obj.loadmcmcdraws(op.Draws);

%% Analyze sample
obj.Post.MCMCSample.Param = sumstats(draws.Param,op.Percentiles);
obj.Post.MCMCSample.AuxParam = sumstats(draws.AuxParam,op.Percentiles);
Var = draws.Param-repmat(obj.Post.Mean,1,draws.N);
Var = Var*Var'/(draws.N-1);
CorrMat = zeros(np);
for jr=1:np
    if Var(jr,jr)==0, continue, end
    for jc=1:np
        if Var(jc,jc)==0, continue, end
        CorrMat(jr,jc) = Var(jr,jc)/sqrt(Var(jr,jr)*Var(jc,jc));
    end
end
obj.Post.Mean = obj.Post.MCMCSample.Param.Mean;
obj.Post.Median = obj.Post.MCMCSample.Param.Median;
obj.Post.SD = obj.Post.MCMCSample.Param.SD;
obj.Post.Prc05 = obj.Post.MCMCSample.Param.Prc05;
obj.Post.Prc95 = obj.Post.MCMCSample.Param.Prc95;
obj.Post.Var = Var;
obj.Post.Corr = CorrMat;


%% check for new mode in mcmc sample
[lpdfNewMode,idxMax] = max(draws.LPDF);
fprintf('Checking MCMC draws for new mode\n')
fprintf('Previous mode lpdf: %.6f\n',obj.Post.ModeLPDF)
fprintf('Highest posterior density in MCMC draws: %.6f\n',lpdfNewMode)
if lpdfNewMode>obj.Post.ModeLPDF
    fprintf('Found MCMC draw with higher posterior density.\n')
    fprintf('Posterior mode updated.\n')
    namelength = [cellfun('length',pNames)];
    namelengthmax = max(namelength);
    fprintf(['%',int2str(namelengthmax),'s %-7s %-7s\n'],'','Old','New')
    for jp=1:np
        fprintf(['%',int2str(namelengthmax),'s %7.4f %7.4f\n'],...
                pNames{jp},obj.Post.Mode(jp),draws.Param(jp,idxMax))
    end
    obj.Post.ModeLPDF = lpdfNewMode;
    obj.Post.Mode = draws.Param(:,idxMax);
else
    fprintf('Did not find MCMC draw with higher posterior density.\n')
    fprintf('Previous posterior mode kept.\n')
end


%% Marginal likelihood
fprintf('Computing marginal likelihood...\n')
tau = 0.1:0.1:0.9;
ntau = length(tau);
pIdx = obj.Post.EstimateIdx;
npd = obj.Post.NEstimate;
% InvVar = inv(Post.Var);
xdd = draws.Param(pIdx,:)-repmat(obj.Post.Mean(pIdx),1,draws.N);
xddvar = xdd*xdd'/draws.N;
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
lpdfMax = obj.Post.ModeLPDF;
lpdfMean = mean(draws.LPDF);
% [postMax,postMean]
% Constant terms
% lfConst = -log(tau)-np/2*log(2*pi)-1/2*log(det(Post.Var));
lfConst = -log(tau)-npd/2*log(2*pi)-1/2*xddvarlndet;
chi2Crit = chi2inv(tau,npd);
% Calculate the ratio of f(x)/post(x)
pw = zeros(ntau,draws.N);
nB = 10;
for jB=1:nB
%     fprintf('Computing %.0f of %.0f...\n',jB,nB)
    for jd=1:draws.N/nB
        idx = ceil(draws.N/nB*(jB-1)+jd);
        pratio = xdd(:,idx)'*InvVar*xdd(:,idx);
        pw(:,idx) = (pratio<=chi2Crit).*...
            exp(lfConst-1/2*pratio-draws.LPDF(idx)+lpdfMean);
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
obj.Post.LogMgLikelihood = LogMgLikelihood(tau==0.5);
% obj.LogMgLikelihood = mean(LogMgLikelihood);
for jt = 1:ntau
    obj.Post.MCMCSample.LogMgLikelihood(jt).Tau = tau(jt);
    obj.Post.MCMCSample.LogMgLikelihood(jt).Value = LogMgLikelihood(jt);
end


%% display results on screen
pNames = obj.Param.Names;
pNameLength = [cellfun('length',pNames)];
nameLengthMax = max(pNameLength);
auxN = obj.AuxParam.N;
auxNames = obj.AuxParam.Names;
auxNameLength = [cellfun('length',auxNames)];
nameLengthMax = max(nameLengthMax,max(auxNameLength));

fprintf('\nResults from MCMC sample analysis:')
fprintf('\n==================================\n')
DispList = {'','',pNames;
            'Prior',' Dist',obj.Prior.Dist;
%             '','   Mode',obj.Prior.Mode;
            '','   Mean',obj.Prior.Mean;
            '','     SD',obj.Prior.SD;
            '','     5%',obj.Prior.Prc05;
            '',' Median',obj.Prior.Median;
            '','    95%',obj.Prior.Prc95;
            'Posterior','   Mode',obj.Post.Mode;
            '','   Mean',obj.Post.Mean;
            '','     SD',obj.Post.SD;
            '','     5%',obj.Post.Prc05;
            '',' Median',obj.Post.Median;
            '','    95%',obj.Post.Prc95;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
    str2show = sprintf('%s  %-5s',str2show,DispList{2,jr});
    for jc=3:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:obj.Param.N
    str2show = sprintf(['%',int2str(nameLengthMax),'s'],DispList{1,3}{jp});
    str2show = sprintf('%s  %5s',str2show,DispList{2,3}{jp});
    for jc=3:nc
        str2show = sprintf('%s  %7.3f',str2show,DispList{jc,3}(jp));
    end
    disp(str2show)
end
fprintf('\n')

DispList = {'','',auxNames;
            'Prior','   Mean',obj.Prior.Sample.AuxParam.Mean;
            '','     SD',obj.Prior.Sample.AuxParam.SD;
            '','     5%',obj.Prior.Sample.AuxParam.Prc05;
            '',' Median',obj.Prior.Sample.AuxParam.Median;
            '','    95%',obj.Prior.Sample.AuxParam.Prc95;
            'Posterior','   Mean',obj.Post.MCMCSample.AuxParam.Mean;
            '','     SD',obj.Post.MCMCSample.AuxParam.SD;
            '','     5%',obj.Post.MCMCSample.AuxParam.Prc05;
            '',' Median',obj.Post.MCMCSample.AuxParam.Median;
            '','    95%',obj.Post.MCMCSample.AuxParam.Prc95;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(nameLengthMax),'s'],DispList{1,jr});
    for jc=2:nc
        str2show = sprintf('%s  %-7s',str2show,DispList{jc,jr});
    end
    disp(str2show)
end
for jp=1:auxN
    str2show = sprintf(['%',int2str(nameLengthMax),'s'],...
                       DispList{1,3}{jp});
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

%% show correlation matrix
if op.showcorr
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
                               obj.Post.Corr(jr,jc));
        end
        disp(str2show)
    end
    disp(' ')
end

%% Make Prior Post Plots
fprintf('Making Prior-Posterior Plots...\n')
fn = sprintf('%s%s-plots-mcmc-%.0f-priorpost',...
             op.PlotDir,obj.Name,obj.Post.MCMCStage);
nPlots = prod([op.Fig.Shape{:}]);
drawsPrior.Param = obj.priordraw(op.NDrawsPrior);
dAux = nan(nAux,op.NDrawsPrior);
parfor jd=1:op.NDrawsPrior
    Matsj = obj.mats(drawsPrior.Param(:,jd));
    dAux(:,jd) = Matsj.AuxParam;
end
drawsPrior.AuxParam = dAux;
drawsPrior.N = size(dAux,2);
pList = {'Param','AuxParam'};
pListPretty = {'Parameters','Auxiliary Parameters'};
dList = {'Prior','Post'};
nFig = ones(1,2);
for jP=1:2
    Pj = pList{jP};
    xj.Post = draws.(Pj);
    xj.Prior = drawsPrior.(Pj);
    xBounds = [min(obj.Prior.Sample.(Pj).Prc01,obj.Post.MCMCSample.(Pj).Prc01),...
               max(obj.Prior.Sample.(Pj).Prc99,obj.Post.MCMCSample.(Pj).Prc99)];
    np = obj.(Pj).N;
    nFig(jP) = ceil(np/nPlots);
    for jF=1:nFig(jP)
        hfig = figure('Visible',op.Fig.Visible);
        clear hf
        for jf=1:nPlots 
            jp = (jF-1)*nPlots+jf;
            if jp>np, break, end
            hf(jf) = subplot(op.Fig.Shape{:},jf);
            for jD=1:2
                Dj = dList{jD};
                xjdata = xj.(Dj)(jp,:);
                xjdata(xjdata>xBounds(jp,2)) = [];
                xjdata(xjdata<xBounds(jp,1)) = [];
                binWidth = xBounds(jp,:)*[-1;1]/op.NBin;
                if binWidth~=0
                    histogram(xjdata,xBounds(jp,1):binWidth:xBounds(jp,2),...
                              'Normalization','probability',...
                              'FaceColor',op.Fig.Color(jD,:),...
                              'FaceAlpha',0.6);
                else
                    histogram(xjdata,11,...
                          'Normalization','probability',...
                          'FaceColor',op.Fig.Color(jD,:),...
                          'FaceAlpha',0.6);
                end
                hold on
            end
            hold off
            axis tight
            ax = gca;
            ax.YTick = [];
            ax.FontSize = op.Fig.FontSize;
            if jf==nPlots || jp==np
                hl = legend(dList,'Orientation','horizontal');
                legPos = hl.Position;
                legPos(1) = 0.5 - legPos(3)/2;
                legPos(2) = 0;
                hl.Position = legPos;
                dLeg = legPos(4);
            end
            title(obj.(Pj).PrettyNames{jp})
        end
        hfig.PaperSize = op.PaperSize;
        hfig.PaperPosition = op.PaperPosition;
        if op.TightFig
            tightfig(hfig,op.Fig.Shape,hf,op.TightFigOptions)
        end
        print('-dpdf',sprintf('%s-%s-%.0f',fn,lower(Pj),jF))
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
fprintf(fid,'thinning used: & %.0f\\\\\n',sample.NThinning);
fprintf(fid,'number of draws used: & %.0f\\\\\\\\\n',draws.N);
fprintf(fid,'log-marginal likelihood: & %.4f\n',obj.Post.LogMgLikelihood);
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{equation*}\n');
fprintf(fid,'\\newpage \n');

fprintf(fid,'\\section{Tables}\n');
% fprintf(fid,'\\subsection{Parameters}\n');
np = obj.Param.N;
str = [' & $%.',int2str(op.Table.Precision),'f$'];
tableBreaks = settablebreaks(np,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\subsection{Parameters}\n');
    else
        fprintf(fid,'\\subsection{Parameters (%.0f/%.0f)}\n',jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    if op.Table.MoveLeft
        fprintf(fid,'\\hspace{-0.5in}\n');
    end
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+4+1+6));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{4}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{6}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& Dist & 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mode & Mean & SD & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.Param.PrettyNames{jr});
        fprintf(fid,' & %s', obj.Prior.Dist{jr});
        fprintf(fid,str,obj.Prior.Prc05(jr));
        fprintf(fid,str,obj.Prior.Median(jr));
        fprintf(fid,str,obj.Prior.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Post.Mode(jr));
        fprintf(fid,str,obj.Post.Mean(jr));
        fprintf(fid,str,obj.Post.SD(jr));
        fprintf(fid,str,obj.Post.Prc05(jr));
        fprintf(fid,str,obj.Post.Median(jr));
        fprintf(fid,str,obj.Post.Prc95(jr));
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

tableBreaks = settablebreaks(auxN,op.Table.MaxRows);
idxPar = 0;
nBreaks = length(tableBreaks);
for jBreak=1:nBreaks
    idxPar = (idxPar(end)+1):tableBreaks(jBreak);
    if nBreaks==1
        fprintf(fid,'\\subsection{Auxiliary Parameters}\n');
    else
        fprintf(fid,'\\subsection{Auxiliary Parameters (%.0f/%.0f)}\n',...
                jBreak,nBreaks);
    end
    fprintf(fid,'\\begin{equation*}\n');
    fprintf(fid,'\\begin{tabular}{l%s} \n',repmat('r',1,1+3+1+5));
    fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
    fprintf(fid,'& \\multicolumn{3}{c}{Prior} ');
    fprintf(fid,'& & \\multicolumn{5}{c}{Posterior} \\\\[0.5ex]\n');
    fprintf(fid,'& 5\\%% & Median & 95\\%% ');
    fprintf(fid,'& & Mean & SD & 5\\%% & Median & 95\\%% \n');
    fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
    for jr=idxPar
        fprintf(fid,'%s',obj.AuxParam.PrettyNames{jr});
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Prior.Sample.AuxParam.Prc95(jr));
        fprintf(fid,' &');
        fprintf(fid,str,obj.Post.MCMCSample.AuxParam.Mean(jr));
        fprintf(fid,str,obj.Post.MCMCSample.AuxParam.SD(jr));
        fprintf(fid,str,obj.Post.MCMCSample.AuxParam.Prc05(jr));
        fprintf(fid,str,obj.Post.MCMCSample.AuxParam.Median(jr));
        fprintf(fid,str,obj.Post.MCMCSample.AuxParam.Prc95(jr));
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

fprintf(fid,'\\section{Prior-Post PDF}\n');
for jP=1:2
    for jF=1:nFig(jP)
        if nFig(jP)==1
            fprintf(fid,'\\subsection{%s}\n',pListPretty{jP});
        else
            fprintf(fid,'\\subsection{%s (%.0f/%.0f)}\n',pListPretty{jP},...
                jF,nFig(jP));
        end
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{Fig_%s_%.0f}\n',pList{jP},jF);
        fprintf(fid,'\\includegraphics[width=\\textwidth]{%s-%s-%.0f.pdf}\n',...
                fn,lower(pList{jP}),jF);
%         fprintf(fid,['\\includegraphics[width=\\textwidth,clip,viewport=' ...
%                      '130 230 490 540]{%s-%s-%.0f.pdf}\n'],fn,pList{jP},jF);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\clearpage \n');
    end
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)





end

