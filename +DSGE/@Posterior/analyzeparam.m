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

op = updateoptions(op,varargin{:});

%% Preparations

fprintf('\n*** Analyzing MCMC Sample %.0f\n',obj.MCMCStage)
ttName = sprintf('AnalyzeParamMCMC%.0f',obj.MCMCStage);
obj.TimeElapsed.start(ttName)

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
xdAux = zeros(nAux,nDrawsUsed);
fh = @(x)obj.Model.mats(x);
parfor jd=1:nDrawsUsed
    Matsj = fh(xd(:,jd));
    xdAux(:,jd) = Matsj.AuxParam;
end


%% Analyze sample
obj.MCMCSample.Param = sumstats(xd,op.Percentiles);
obj.MCMCSample.AuxParam = sumstats(xdAux,op.Percentiles);
obj.Mean = obj.MCMCSample.Param.Mean;
obj.Median = obj.MCMCSample.Param.Median;
obj.SD = obj.MCMCSample.Param.SD;
obj.Prc05 = obj.MCMCSample.Param.Prc05;
obj.Prc95 = obj.MCMCSample.Param.Prc95;
Var = xd-repmat(obj.Mean,1,nDrawsUsed);
obj.Var = Var*Var'/(nDrawsUsed-1);
CorrMat = zeros(np);
for jr=1:np
    for jc=1:np
        CorrMat(jr,jc) = obj.Var(jr,jc)/sqrt(obj.Var(jr,jr)*obj.Var(jc,jc));
    end
end
obj.Corr = CorrMat;

%% check for new mode in mcmc sample
[lpdfNewMode,idxMax] = max(lpdfd);
fprintf('\nChecking MCMC draws for new mode\n')
fprintf('Previous mode lpdf: %.6f\n',obj.LPDFMode)
fprintf('Highest posterior density in MCMC draws: %.6f\n',lpdfNewMode)
if lpdfNewMode>obj.LPDFMode
    fprintf('Found MCMC draw with higher posterior density.\n')
    fprintf('Posterior mode updated!\n')
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
    fprintf('Previous posterior mode kept!\n')
end
fprintf('\n')

%% Marginal likelihood
fprintf('Computing marginal likelihood\n')
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
            'Prior','Dist',obj.Prior.Dist;
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
            'Posterior','   Mode',xAux;
            '','   Mean',obj.Sample.AuxParam.Mean;
            '','     SD',obj.Sample.AuxParam.SD;
            '','     5%',obj.Sample.AuxParam.Prc05;
            '',' Median',obj.Sample.AuxParam.Median;
            '','    95%',obj.Sample.AuxParam.Prc95;
           };
nc = size(DispList,1);
for jr=1:2
    str2show = sprintf(['%-',int2str(pNameLengthMax),'s'],DispList{1,jr});
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

%% Finish up
obj.TimeElapsed.stop(ttName)

