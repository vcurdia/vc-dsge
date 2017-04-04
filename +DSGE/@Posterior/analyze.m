function analyze(obj,varargin)

% analyze
% 
% Analyze posterior MCMC sample 
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.SampleID = length(obj.Sample);
op.BurnIn = 0.5;
op.Thinning = 1;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op.Table = DSGE.Options.Table;

op = updateoptions(op,varargin{:});

%% Preparations

sid = op.SampleID;
sample = obj.Sample(sid);
pN = obj.Model.Param.N;
pNames = obj.Model.Param.Names;
auxN = obj.Model.AuxParam.N;

fprintf('\n*** Analyzing MCMC Sample %.0f\n',sid)
ttName = sprintf('AnalyzeSample%.0f',sid);
obj.TimeElapsed.start(ttName)

%% load the mcmc draws
idxDraws = (op.BurnIn*sample.NDraws+1):op.NThinning:sample.NDraws;
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
fprintf('Thinning used: %.0f\n', op.NThinning)
fprintf('Total number of draws used: %.0f\n', nDrawsUsed)

%% Generate AuxParam
xdAux = zeros(auxN,nDrawsUsed);
fh = @(x)obj.Model.mats(x);
parfor jd=1:nDrawsUsed
    Matsj = fh(xd(:,jd));
    xdAux(:,jd) = Matsj.AuxParam;
end


%% Analyze sample
obj.Sample(sid).Param = sumstats(xd,op.Percentiles);
obj.Sample(sid).AuxParam = sumstats(xdAux,op.Percentiles);
obj.Mean = obj.Sample(sid).Param.Mean;
obj.Median = obj.Sample(sid).Param.Median;
obj.SD = obj.Sample(sid).Param.SD;
obj.Prc05 = obj.Sample(sid).Param.Prc05;
obj.Prc95 = obj.Sample(sid).Param.Prc95;
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

Mg likelihood next


%% Finish up
obj.TimeElapsed.stop(ttName)

