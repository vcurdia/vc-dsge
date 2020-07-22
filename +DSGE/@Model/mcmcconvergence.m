function mcmcconvergence(obj,varargin)

% mcmcconvergence
% 
% Analyze convergence of MCMC sample
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 8, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% Options
op.Draws.BurnIn = 0.25;
op.Draws.AuxParam = 0;
op.NMeansSPM = 4;
op.DrawsFraction = 0.04;
op.Table = DSGE.Options.Table;
op.NBin = 50;
op.TraceStep = [];
op.Fig.Visible = 'off';
op.Fig.Color = colorscheme;
op.Fig.FontSize = 6;
op.PaperSize = [6.5, 6.5];
op.PaperPosition = [0. 0, 6.5, 6.5];
op.TightFig = 1;
op.TightFigOptions = struct;
op.PlotDirDraws = 'plots-draws/';
op.PlotDirTrace = 'plots-trace/';

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\nAnalyzing convergence of MCMC Sample %.0f\n',obj.Post.MCMCStage)

if ~isdir(op.PlotDirDraws),mkdir(op.PlotDirDraws),end
if ~isdir(op.PlotDirTrace),mkdir(op.PlotDirTrace),end
ReportFileName = sprintf('report-%s-mcmc-conv-%.0f',obj.Name,...
                         obj.Post.MCMCStage);
ReportTitle = sprintf('%s\\\\[30pt]Convergence Analysis\\\\MCMC Stage %.0f',...
                      obj.Name,obj.Post.MCMCStage);

sample = obj.Post.MCMCSample;

%% load the mcmc draws
draws = obj.loadmcmcdraws(op.Draws,'BurnIn',0,'CombineChains',0,'ExpandParam',0);
nDrawsUsed = draws.N;
idxburned = ceil(op.Draws.BurnIn*nDrawsUsed+1:nDrawsUsed);


if isempty(op.NBin), op.NBin = round(2*nDrawsUsed^(1/3)); end

p.LPDF.Title = 'Log-PDF';
p.LPDF.Names = {'LPDF'};
p.LPDF.PrettyNames = {p.LPDF.Title};
p.Param.Title = 'Parameters';
p.Param.Names = obj.Param.Names(obj.Post.EstimateIdx);
p.Param.PrettyNames = obj.Param.PrettyNames(obj.Post.EstimateIdx);
p.AuxParam.Title = 'Auxiliary Parameters';
p.AuxParam.Names = obj.AuxParam.Names;
p.AuxParam.PrettyNames = obj.AuxParam.PrettyNames;
pList = fieldnames(p);
nList = length(pList);
for jL=1:nList
    Lj = pList{jL};
    p.(Lj).N = length(p.(Lj).Names);
    p.(Lj).NameLength = [cellfun('length',p.(Lj).Names)];
    p.(Lj).NameLengthMax = max(p.(Lj).NameLength);
end
nameLengthMax = p.Param.NameLengthMax;
if op.Draws.AuxParam
    nameLengthMax = max(nameLengthMax,p.AuxParam.NameLengthMax);
end


%% Plot draws
fprintf('Making plots of MCMC draws\n')
pdFN = sprintf('%s%s-plots-mcmc-%.0f-draws',...
             op.PlotDirDraws,obj.Name,obj.Post.MCMCStage);
if op.Draws.AuxParam
    pdList = {'LPDF','Param','AuxParam'};
else
    pdList = {'LPDF','Param'};
end    
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
        hf = figure('Visible',op.Fig.Visible);
        clear h
        jh = 0;
        for jChain=1:sample.NChains
            jh = jh+1;
            h(jh) = subplot(sample.NChains,nc,jh);
            plot(xd(1,:,jChain),'LineWidth',1)
            if pMax>pMin
                ylim([pMin pMax])
            end
            h(jh).XTick = [op.Draws.BurnIn,1]*nDrawsUsed;
            h(jh).XGrid = 'on';
            h(jh).FontSize = op.Fig.FontSize;
            if nc>1
                jh = jh+1;
                h(jh) = subplot(sample.NChains,2,jh);
                histogram(xd(1,idxburned,jChain),...
                          op.NBin,'Normalization','probability',...
                          'FaceColor',op.Fig.Color(1,:),'FaceAlpha',1)
                if pMax>pMin
                    xlim([pMin pMax])
                end
                h(jh).FontSize = op.Fig.FontSize;
            end
        end
        title(h(1),sprintf('%s in each chain',p.(Lj).PrettyNames{jF}))
        if nc>1
            title(h(2),sprintf('Hist excluding initial %.0f\\%% of obs',...
                                 100*op.Draws.BurnIn))
        end
        hf.PaperSize = op.PaperSize;
        hf.PaperPosition = op.PaperPosition;
        if op.TightFig
            tightfig(hf,{sample.NChains,nc},h,op.TightFigOptions)
        end
        print('-dpdf',sprintf('%s-%s.pdf',pdFN,p.(Lj).Names{jF}))
    end
end
close all

%% eliminate BurnIn
draws.Param = draws.Param(:,idxburned,:);
draws.LPDF = draws.LPDF(:,idxburned,:);
if op.Draws.AuxParam
    draws.AuxParam = draws.AuxParam(:,idxburned,:);
end
nDrawsUsed = size(draws.LPDF,2);
fprintf('Burn in for rest of convergence analysis: %.0f%%\n',...
        100*op.Draws.BurnIn)

%% Plot trace
fprintf('Making trace plots\n')
ptFN = sprintf('%s%s-plots-mcmc-%.0f-trace',...
             op.PlotDirTrace,obj.Name,obj.Post.MCMCStage);
if op.Draws.AuxParam
    ptList = {'Param','AuxParam'};
else
    ptList = {'Param'};
end    
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
        hf = figure('Visible',op.Fig.Visible);
        clear h
        jh=0;
        for jChain=1:sample.NChains
            jh = jh+1;
            h(jh) = subplot(sample.NChains,2,jh);
            plot(SampleID,pMean*ones(size(SampleID)),'-',...
                 'Color',op.Fig.Color(1,:))
            hold on
            plot(SampleID,(pMean-pSD)*ones(size(SampleID)),':',...
                 'Color',op.Fig.Color(2,:))
            plot(SampleID,(pMean+pSD)*ones(size(SampleID)),':',...
                 'Color',op.Fig.Color(2,:))
            plot(SampleID,RollingMean(:,jChain),'-','Color',op.Fig.Color(1,:),...
                 'LineWidth',2)
            if MeanBounds(2)>MeanBounds(1)
                ylim(MeanBounds)
            end
            xlim(SampleID([1,end]))
            h(jh).FontSize = op.Fig.FontSize;
            jh = jh+1;
            h(jh) = subplot(sample.NChains,2,jh);
            plot(SampleID,pSD*ones(size(SampleID)),'-',...
                 'Color',op.Fig.Color(1,:))
            hold on
            plot(SampleID,RollingSD(:,jChain),'-','Color',op.Fig.Color(1,:),...
                 'LineWidth',2)
            if SDBounds(2)>SDBounds(1)
                ylim(SDBounds)
            end
            xlim(SampleID([1,end]))
            h(jh).FontSize = op.Fig.FontSize;
        end
        title(h(1),sprintf('Rolling Mean of %s',p.(Lj).PrettyNames{jF}))
        title(h(2),sprintf('Rolling SD of %s',p.(Lj).PrettyNames{jF}))
        hf.PaperSize = op.PaperSize;
        hf.PaperPosition = op.PaperPosition;
        if op.TightFig
            tightfig(hf,{sample.NChains,2},h,op.TightFigOptions)
        end
        print('-dpdf',sprintf('%s-%s.pdf',ptFN,p.(Lj).Names{jF}))
    end
end
close all

%% Convergence Tests
if op.Draws.AuxParam
    cList = {'Param','AuxParam'};
else
    cList = {'Param'};
end

%% R and MNEff from Gelman et al
fprintf('\nR statistic from Gelman et al:')
fprintf('\n==============================\n\n')
for jL=1:length(cList)
    Lj = cList{jL};
    xd = draws.(Lj);
    Mj = mean(xd,2);
    M = mean(Mj,3);
    s2j = sum((xd-repmat(Mj,[1,nDrawsUsed,1])).^2,2)/(nDrawsUsed-1);
    W = mean(s2j,3);
    B = nDrawsUsed/(sample.NChains-1)*...
        sum((Mj-repmat(M,[1,1,sample.NChains])).^2,3);
    varplus = (nDrawsUsed-1)/nDrawsUsed*W+1/nDrawsUsed*B;
    conv.R.(Lj) = squeeze(sqrt(varplus./W));
    conv.MNEff.(Lj) = squeeze(sample.NChains*nDrawsUsed*varplus./B);
    for jp=1:p.(Lj).N
        fprintf(['%',int2str(nameLengthMax),'s %8.4f\n'],...
                p.(Lj).Names{jp},conv.R.(Lj)(jp))
    end
    fprintf('\n')
end
fprintf('\nEffective number of independent draws a la Gelman et al.:')
fprintf('\n=========================================================\n\n')
for jL=1:length(cList)
    Lj = cList{jL};
    for jp=1:p.(Lj).N
        fprintf(['%',int2str(nameLengthMax),'s %10.0f\n'],...
                 p.(Lj).Names{jp},conv.MNEff.(Lj)(jp))
    end
    fprintf('\n')
end


%% Geweke's separated partial means test
npm = floor(nDrawsUsed/2/op.NMeansSPM);
Lp = round(op.DrawsFraction*npm);
SPMc95 = chi2inv(0.95,op.NMeansSPM-1);
SPMc99 = chi2inv(0.99,op.NMeansSPM-1);
fprintf('\nSPM test results for each chain:')
fprintf('\n================================')
fprintf('\ncritical value (95%%) for chi-square(%.0f) is %f\n\n',...
        op.NMeansSPM-1,chi2inv(0.95,op.NMeansSPM-1))
for jL=1:length(cList)
    Lj = cList{jL};
    xd = draws.(Lj);
    SPMj = zeros(p.(Lj).N,sample.NChains);
    parfor jChain=1:sample.NChains
        mean_jp = zeros(p.(Lj).N,op.NMeansSPM);
        S0 = zeros(p.(Lj).N,op.NMeansSPM);
        for jp=1:op.NMeansSPM
            xj = xd(:,(2*jp-1)*npm+1:2*jp*npm,jChain);
            mean_jp(:,jp) = mean(xj,2);
            xj = xj-repmat(mean_jp(:,jp),1,npm);
            s0 = mean(xj.^2,2);
            for jL=1:Lp-1
                cL = sum(xj(:,1+jL:npm).*xj(:,1:npm-jL),2)/npm;
                s0 = s0 + 2*(Lp-jL)/Lp*cL;
            end
            S0(:,jp) = s0; 
        end
        SPMjj = zeros(p.(Lj).N,1);
        for j=1:p.(Lj).N
            hp = mean_jp(j,2:end)'-mean_jp(j,1:end-1)';
            Vp = diag(S0(j,2:end))+diag(S0(j,1:end-1));
            Vp = Vp - [zeros(op.NMeansSPM-1,1),...
                       [diag(S0(j,2:end-1));zeros(1,op.NMeansSPM-2)]];
            Vp = Vp - [zeros(1,op.NMeansSPM-1);...
                       [diag(S0(j,2:end-1)),zeros(op.NMeansSPM-2,1)]];
            Vp = Vp/npm;
            SPMjj(j) = hp'*inv(Vp)*hp;
        end
        SPMj(:,jChain) = SPMjj;
    end
    conv.SPM.(Lj) = SPMj;
    for jp=1:p.(Lj).N
        fprintf(['%',int2str(nameLengthMax),'s', ...
                 repmat(' %10.4f',1,sample.NChains),'\n'],...
                p.(Lj).Names{jp},conv.SPM.(Lj)(jp,:))
    end
    fprintf('\n')
end

%% Calculate the within chain number of effective sample size
L = round(op.DrawsFraction*nDrawsUsed);
fprintf('\nn_eff for each chain:')
fprintf('\n=====================\n\n')
for jL=1:length(cList)
    Lj = cList{jL};
    xd = draws.(Lj);
    neffj = zeros(p.(Lj).N,sample.NChains);
    parfor jChain=1:sample.NChains
        xj = xd(:,:,jChain);
        xj = xj-repmat(mean(xj,2),1,nDrawsUsed);
        c0 = mean(xj.^2,2);
        S0 = c0;
        for jL=1:L-1
            cL = sum(xj(:,1+jL:nDrawsUsed).*xj(:,1:nDrawsUsed-jL),2)/nDrawsUsed;
            S0 = S0 + 2*(L-jL)/L*cL;
        end
        neffj(:,jChain) = nDrawsUsed*c0./S0;
    end
    conv.NEff.(Lj) = neffj;
    for jp=1:p.(Lj).N
        fprintf(['%',int2str(nameLengthMax),'s',...
                 repmat(' %10.0f',1,sample.NChains),'\n'],...
                p.(Lj).Names{jp},conv.NEff.(Lj)(jp,:))
    end
    fprintf('\n')
end

obj.Post.MCMCSample.Convergence = conv;

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
fprintf(fid,'\\newpage\n');

fprintf(fid,'\\section{R, Effective Draws}\n');
for jL=1:length(cList)
    Lj = cList{jL};
    np = p.(Lj).N;
    tableBreaks = settablebreaks(np,op.Table.MaxRows);
    idxPar = 0;
    nBreaks = length(tableBreaks);
    for jBreak=1:nBreaks
        idxPar = (idxPar(end)+1):tableBreaks(jBreak);
        if nBreaks==1
            fprintf(fid,'\\subsection{%s}\n',p.(Lj).Title);
        else
            fprintf(fid,'\\subsection{%s (%.0f/%.0f)}\n',p.(Lj).Title,...
                    jBreak,nBreaks);
        end
        fprintf(fid,'\\begin{equation*}\n');
        if op.Table.MoveLeft
            fprintf(fid,'\\hspace{-0.5in}\n');
        end
        fprintf(fid,'\\begin{tabular}{lrr%s} \n',repmat('r',1,sample.NChains));
        fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
        fprintf(fid,['& & \\multicolumn{%.0f}{c}{Number of Effective ' ...
                     'Draws}\\\\[0.5ex]\\cline{3-%.0f}\\\\[-1.5ex]\n'],...
                1+sample.NChains,3+sample.NChains);
        fprintf(fid,'& \\multicolumn{1}{c}{$R$} & Sample');
        for jChain=1:sample.NChains
            fprintf(fid,' & Chain $%.0f$\n',jChain);
        end
        fprintf(fid,'\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
        for jr=idxPar
            fprintf(fid,'%s',p.(Lj).PrettyNames{jr});
            fprintf(fid,' & $%.4f$',conv.R.(Lj)(jr));
            fprintf(fid,' & $%.0f$',conv.MNEff.(Lj)(jr));
            for jChain=1:sample.NChains
                fprintf(fid,' & $%.0f$',conv.NEff.(Lj)(jr,jChain));
            end
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
end

fprintf(fid,'\\section{SPM}\n');
for jL=1:length(cList)
    Lj = cList{jL};
    np = p.(Lj).N;
    tableBreaks = settablebreaks(np,op.Table.MaxRows);
    idxPar = 0;
    nBreaks = length(tableBreaks);
    for jBreak=1:nBreaks
        idxPar = (idxPar(end)+1):tableBreaks(jBreak);
        if nBreaks==1
            fprintf(fid,'\\subsection{%s}\n',p.(Lj).Title);
        else
            fprintf(fid,'\\subsection{%s (%.0f/%.0f)}\n',p.(Lj).Title,...
                    jBreak,nBreaks);
        end
        fprintf(fid,'\\begin{equation*}\n');
        if op.Table.MoveLeft
            fprintf(fid,'\\hspace{-0.5in}\n');
        end
        fprintf(fid,'\\begin{tabular}{l%s}\n',repmat('rl',1,sample.NChains));
        fprintf(fid,'\\hline\\hline\\\\[-1.5ex]\n');
        for jChain=1:sample.NChains
            fprintf(fid,' & \\multicolumn{2}{c}{$SPM_%.0f(%.0f)$}',...
                    op.NMeansSPM,jChain);
        end
        fprintf(fid,'\n\\\\[0.5ex]\\hline\\\\[-1.5ex]\n');
        for jr=idxPar
            fprintf(fid,'%s',p.(Lj).PrettyNames{jr});
            for jChain=1:sample.NChains
                fprintf(fid,' & $%.2f$',conv.SPM.(Lj)(jr,jChain));
                fprintf(fid,' & %s%s',...
                        repmat('*',conv.SPM.(Lj)(jr,jChain)>=SPMc99),...
                        repmat('*',conv.SPM.(Lj)(jr,jChain)>=SPMc95));
            end
            fprintf(fid,' \\\\\n');
            if ismember(jr,op.Table.Lines) && jr~=idxPar(end)
                fprintf(fid,'\\\\[-1.5ex]\\hline\\\\[-1.5ex]\n');
            end        
        end
        fprintf(fid,'\\\\[-1.5ex]\\hline\\hline\n');
        fprintf(fid,'\\end{tabular}\n');
        fprintf(fid,'\\end{equation*}\n');
        fprintf(fid,'The null hypothesis of the SPM test is that the mean \n');
        fprintf(fid,'in two separate subsamples is the same. * indicates \n');
        fprintf(fid,'p-value less than 5\\%%. ** indicates p-value less \n');
        fprintf(fid,'than 1\\%%.\n');
        fprintf(fid,'\\clearpage\n');
    end
end

fprintf(fid,'\\section{Draws}\n');
for jL=1:length(pdList)
    Lj = pdList{jL};
    fprintf(fid,'\\subsection{%s}\n',p.(Lj).Title);
    for jF=1:p.(Lj).N
        if p.(Lj).N>1
            fprintf(fid,'\\subsubsection{%s}\n',p.(Lj).Names{jF});
        end
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\label{Fig_%s}\n',p.(Lj).Names{jF});
        fprintf(fid,'\\includegraphics[width=\\textwidth]{%s-%s.pdf}\n',...
                pdFN,p.(Lj).Names{jF});
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
        fprintf(fid,'\\includegraphics[width=\\textwidth]{%s-%s.pdf}\n',...
                ptFN,p.(Lj).Names{jF});
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\clearpage \n');
    end
end

fprintf(fid,'\\end{document}\n');
fclose(fid);
pdflatex(ReportFileName)



end

