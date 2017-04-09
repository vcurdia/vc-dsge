function nRejections = mcmcchain(obj,varargin)

% mcmcchain
% 
% Generate MCMC chain
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 28, 2017
% Copyright (C) 2017 Vasco Curdia

%% map object related variables
np = obj.NEstimate;
xMode = obj.Mode(obj.EstimateIdx);
Var = obj.Var(obj.EstimateIdx,obj.EstimateIdx);
lpdfMode = obj.LPDFMode;

%% Options
op.verbose = 1;
op.Augment = 0;
op.NDraws = 50000;
op.NIRS = 1000;
op.NBlocks = 10;
op.ExplosionScale = [.4,.5,.6,.7];
op.ExplosionProb = [.3,.6,.8,1];
op.NChangedNIRS = 2;
op.InitDrawTStudent = 1;
op.InitDrawItMax = 1000;
op.InitDrawTol = 20;
op.InitDrawVarFactor = 1;
op.InitDrawDF = 6;
op.JumpVar = 2.4^2/np*Var;
op.fn = 'MCMC_Chain';
op.x0 = [];

op = updateoptions(op,varargin{:});

%% chain related variables
fid = fopen([op.fn,'.log'],'wt');
lpdf = @(x)obj.lpdf(x,struct('verbose',op.verbose,'fid',fid));

%% generate initial draw
x0 = op.x0;
if ~op.Augment && isempty(x0)
    fprintf(fid,'Generating initial draw...\n');
    if op.InitDrawTStudent
        InitDrawVarChol = chol(op.InitDrawVarFactor^2*Var)';
        InitDrawCorr = eye(np);
        for jd=1:op.InitDrawItMax
            x0 = xMode + InitDrawVarChol*mvtrnd(InitDrawCorr,op.InitDrawDF)';
            lpdf0 = lpdf(x0);
            if lpdf0>lpdfMode-op.InitDrawTol
                break
            end
        end
        if lpdf0<lpdfMode-op.InitDrawTol
            fprintf(fid,...
                    'Warning: Initial draw has low posterior density.\n');
        end
    else
        nExplosion = length(op.ExplosionScale);
        nIRS = op.NIRS;
        jChangedNIRS = 0;
        while jChangedNIRS<=op.NChangedNIRS
            ExplosionDraws = unifrnd(0,1,1,nIRS);
            x0c = zeros(np,nIRS);
            x0c1 = zeros(np,nIRS);
            lpdf0c = zeros(1,nIRS);
            lpdf0c1 = zeros(1,nIRS);
            for jd=1:nIRS
                for je=1:nExplosion
                    if ExplosionDraws(jd)<=op.ExplosionProb(je)
                        ExpV = Var*op.ExplosionScale(je)^2;
                        break
                    end
                end
                x0c(:,jd) = mvnrnd(xMode,ExpV)';
                lpdf0c(jd) = lpdf(x0c(:,jd));
            end
            lpdf0c = exp(lpdf0c - repmat(lpdfMode,1,nIRS));
            lpdf0c = lpdf0c./sum(lpdf0c);
            [clpdf0c,idx0] = sort(lpdf0c);
            clpdf0c = cumsum(clpdf0c);
            pick0 = unifrnd(0,1);
            for j=1:nIRS
                if clpdf0c(j)>=pick0
                    x0 = x0c(:,idx0(j));
                    % check whether the draw was ok...
                    WeightPick = lpdf0c(idx0(j));
                    WeightMax = max(lpdf0c);
                    fprintf(fid,'weight of the picked draw: %.6f\n', ...
                            WeightPick);
                    fprintf(fid,'maximum weight for the draw: %.6f\n', ...
                            WeightMax);
                    break
                end
            end
            if ~isempty(x0)
                break
            else
                if jChangedNIRS==op.NChangedNIRS
                    x0 = pMode;
                    WeightPick = NaN;
                    WeightMax = NaN;
                    fprintf(fid,['Warning: Did not find any suitable candidate ' ...
                                 'after increasing nIRS. Using mode.\n']);
                else
                    nIRS = nIRS*10;
                    jChangedNIRS = jChangedNIRS+1;
                    fprintf(fid,['Warning: Did not find any suitable candidate. ' ...
                                 'Trying with increased nIRS.\n']);
                end
            end
        end
    end
end

%% Prepare variables
if ~op.Augment
    draws.N = 0;
    draws.Param = [];
    draws.LPDF = [];
    draws.NRejections = 0;
else
    fprintf(fid,'Loading existing chain...\n');
    load(op.fn)
    nDraws = op.NDraws - draws.N;
    x0 = draws.Param(:,end);
end
nRejections = draws.NRejections;
nDrawsBlock = ceil(nDraws/op.NBlocks);

%% MCMC
lpdf0 = lpdf(x0);
for jB=1:op.NBlocks
    fprintf(fid,'Generating set %2.0f out of %.0f...\n',jB,op.NBlocks);
    nB = min(nDrawsBlock,nDraws-(jB-1)*nDrawsBlock);
    xB = zeros(np,nB);
    lpdfB = zeros(1,nB);
    for j=1:nB
        xc = mvnrnd(x0,op.JumpVar,1)';
        lpdfc = lpdf(xc);
        if unifrnd(0,1)<exp(lpdfc-lpdf0)
            x0 = xc;
            lpdf0 = lpdfc;
        else
            nRejections = nRejections+1;
        end
        xB(:,j) = x0;
        lpdfB(j) = lpdf0;
    end
    draws.N = draws.N+NB;
    draws.Param = [draws.Param,xB];
    draws.LPDF = [draws.LPDF,lpdfB];
    draws.NRejections = nRejections;
    save(op.fn,'draws')
end

%% show number of rejections
fprintf(fid,'%.0f rejections out of %.0f draws (%.2f%%).\n',...
        nRejections,draws.N,nRejections/draws.N*100);

%% save output
save(op.fn,op.SaveList{:});

%% close printed output file
if fid~=1,fclose(fid);end

