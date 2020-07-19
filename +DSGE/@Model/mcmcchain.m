function nRejections = mcmcchain(obj,varargin)

% mcmcchain
% 
% Generate MCMC chain
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: March 28, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% map object related variables
np = obj.Post.NEstimate;
xMode = obj.Post.Mode(obj.Post.EstimateIdx);
Var = obj.Post.Var(obj.Post.EstimateIdx,obj.Post.EstimateIdx);
lpdfMode = obj.Post.ModeLPDF;

%% Options
op.verbose = 1;
op.Augment = 0;
op.NDraws = 50000;
op.NDrawsKeep = 200000;
op.NIRS = 1000;
op.NBlocks = 10;
op.ExplosionScale = [.4,.5,.6,.7];
op.ExplosionProb = [.3,.6,.8,1];
op.NChangedNIRS = 2;
op.InitDrawTStudent = 1;
op.InitDrawItMax = 1000;
op.InitDrawTol = 50;
op.InitDrawVarFactor = 1;
op.InitDrawDF = 6;
op.InitDrawScale = [1,0.5,0.1];
op.JumpVar = 2.4^2/np*Var;
op.fn = 'MCMC_Chain';
op.x0 = [];
op = updateoptions(op,varargin{:});


%% chain related variables
fid = fopen([op.fn,'.log'],'wt');
lpdf = @(x)obj.postlpdf(x,struct('verbose',op.verbose,'fid',fid));

%% generate initial draw
x0 = op.x0;
if ~op.Augment && isempty(x0)
    fprintf(fid,'Generating initial draw...\n');
    if op.InitDrawTStudent
        InitDrawVarChol = chol(op.InitDrawVarFactor^2*Var)';
        InitDrawCorr = eye(np);
        isBad = 1;
        for jScale=1:length(op.InitDrawScale)
            for jd=1:op.InitDrawItMax
                x0 = xMode + op.InitDrawScale(jScale)*InitDrawVarChol*...
                     mvtrnd(InitDrawCorr,op.InitDrawDF)';
                lpdf0 = lpdf(x0);
                isGood = (lpdf0>lpdfMode-op.InitDrawTol);
                if isGood, break, end
            end
            if isGood, break, end
        end
        if ~isGood
            fprintf(fid,'Warning: Initial draw has low posterior density\n');
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
nDrawsKeep = min(op.NDrawsKeep,op.NDraws);
if op.Augment
    fprintf(fid,'Loading existing chain...\n');
    draws = load(op.fn);
    if ~ismember('NKeep',fieldnames(draws))
        draws.NKeep = draws.N;
        draws.NThinning = draws.N/draws.NKeep;
    end
    nDrawsOld = size(draws.Param,2)*draws.NThinning;
    nOldKeep = nDrawsOld/op.NDraws*nDrawsKeep;
    nThinningOld = nDrawsOld/nOldKeep/draws.NThinning;
    idxold = ceil(nThinningOld:nThinningOld:(nDrawsOld/draws.NThinning));
    draws.Param = draws.Param(:,idxold);
    draws.LPDF = draws.LPDF(:,idxold);
    draws.N = nDrawsOld;
    draws.NKeep = nDrawsKeep;
    draws.NThinning = op.NDraws/nDrawsKeep;
    save(op.fn,'-struct','draws');
    nDraws = op.NDraws - nDrawsOld;
    nNewThinning = nDraws/(nDrawsKeep-nOldKeep);
    x0 = draws.Param(:,end);
    lpdf0 = lpdf(x0);
    fprintf(fid,'Adding %.0f draws to existing chain\n',nDraws);
else
    draws.N = 0;
    draws.Param = [];
    draws.LPDF = [];
    draws.NRejections = 0;
    draws.NKeep = nDrawsKeep;
    draws.NThinning = op.NDraws/nDrawsKeep;
    nDraws = op.NDraws;
    nNewThinning = nDraws/nDrawsKeep;
    lpdf0 = lpdf(x0);
    pNames = obj.Param.Names(obj.Post.EstimateIdx);
    fprintf(fid,'Initial draw:\n');
    nameLength = max([cellfun('length',pNames)]);
    for jp=1:np
        fprintf(fid,['%',int2str(nameLength),'s %7.4f\n'],pNames{jp},x0(jp));
    end
    fprintf(fid,'\nInitial posterior level: %.8f\n\n',lpdf0);
end
nRejections = draws.NRejections;

%% MCMC
nDrawsBlock = ceil(nDraws/op.NBlocks);
for jB=1:op.NBlocks
    nB = min(nDrawsBlock,nDraws-(jB-1)*nDrawsBlock);
    idxDraws = ceil(nNewThinning:nNewThinning:nB);
    njkeep = length(idxDraws);
    xB = zeros(np,njkeep);
    lpdfB = zeros(1,njkeep);
    for j=1:nB
        xc = mvnrnd(x0,op.JumpVar,1)';
        lpdfc = lpdf(xc);
        if unifrnd(0,1)<exp(lpdfc-lpdf0)
            x0 = xc;
            lpdf0 = lpdfc;
        else
            nRejections = nRejections+1;
        end
        [tf,idxj] = ismember(j,idxDraws);
        if tf
            xB(:,idxj) = x0;
            lpdfB(idxj) = lpdf0;
        end
    end
    draws.N = draws.N+nB;
    draws.Param = [draws.Param,xB];
    draws.LPDF = [draws.LPDF,lpdfB];
    draws.NRejections = nRejections;
    save(op.fn,'-struct','draws')
    fprintf(fid,'completed %3.0f%%\n',jB/op.NBlocks*100);
end

%% show number of rejections
fprintf(fid,'%.0f rejections out of %.0f draws (%.2f%%).\n',...
        nRejections,draws.N,nRejections/draws.N*100);

% %% save output (not needed because we save at end of each block)
% save(op.fn,'draws');

%% close printed output file
if fid~=1,fclose(fid);end


end

