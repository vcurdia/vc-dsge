function jumpScale = calibratemcmc(obj,varargin)

% calibratemcmc
% 
% Calibrate MCMC jump distribution
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017-2018 Vasco Curdia

%% Options
op.x0 = [];
op.UsePostDraw = 1;
op.NChains = 4;
op.JumpScale = 2.4;
op.JumpNewVarWeight = 1;
op.NDrawsCalibrate = 1000;
op.KeepFilesCalibrate = 0; 
op.ScaleIncrements = 0.1; %[0.2,0.05,0.01];
op.NConfirm = 1;
op.RejectionRateMax = 0.85;
op.RejectionRateMin = 0.65; 
op.MaxReverseDirection = 2; 
op.MinSearchScale = 0.1; 

op = updateoptions(op,varargin{:});

%% Preparations

if isempty(obj.Post.MCMCStage), obj.Post.MCMCStage = 1; end
fprintf('\nCalibrating Jump Distribution for MCMC Sample %.0f\n', ...
        obj.Post.MCMCStage)

pIdx = obj.Post.EstimateIdx;
jumpScale = op.JumpScale;
jumpVarRaw = obj.Post.Var(pIdx,pIdx)/obj.Post.NEstimate;
if op.JumpNewVarWeight<1 && obj.Post.MCMCStage>1
    jumpVarRaw = op.JumpNewVarWeight^2*jumpVarRaw + ...
        (1-op.JumpNewVarWeight)^2*obj.Post.MCMCSample.JumpVar/...
        obj.Post.MCMCSample.JumpScale^2;
end

opChain.Augment = 0;
opChain.NDraws = op.NDrawsCalibrate;

[npx0,nx0] = size(op.x0);
x0 = cell(1,op.NChains);
if nx0>0
    if npx0==obj.Param.N
        op.x0 = op.x0(obj.Post.EstimateIdx,:); 
    end
    for j=1:nx0
        x0{j} = op.x0(:,j);
    end
end
if op.UsePostDraw && obj.Post.MCMCStage>1
    x0d = obj.postdraw(op.NChains-nx0);
    for j=(nx0+1):op.NChains
        x0{j} = x0d(obj.Post.EstimateIdx,j-nx0);
    end
end

fn = cell(op.NChains,1);
for jChain=1:op.NChains
    fn{jChain} = sprintf('%s-mcmc-%.0f-chain-%.0f-calibratejump',...
                         obj.Name,obj.Post.MCMCStage,jChain);
end

%% Calibrate jump scale
jConfirm = 0;
jIncrement = 1;
nRevertDirection = 0;
currChange = 0;
nBlocks = 0;
RejectionRates = zeros(1,op.NChains);
while jConfirm<=op.NConfirm 
    nBlocks = nBlocks+1;
    opChain.JumpVar = jumpScale^2*jumpVarRaw;
    RejectionRates(nBlocks,:) = zeros(1,op.NChains);
    parfor jChain=1:op.NChains
        opj = opChain;
        opj.fn = fn{jChain};
        opj.x0 = x0{jChain};
        nRejections = obj.mcmcchain(opj);
        RejectionRates(nBlocks,jChain) = nRejections/opj.NDraws;
    end
    fprintf('\nLast chain results:\n')
    nAbove = 0;
    nBelow = 0;
    for jChain=1:op.NChains
        JumpScaleFactors(nBlocks) = jumpScale;
        fprintf(['Block %03.0f, Chain %02.0f: JumpScale = %4.2f, ',...
                 'rejection rate = %5.1f%%\n'],...
                nBlocks,jChain,jumpScale,RejectionRates(nBlocks,jChain)*100)
        nAbove = nAbove + (RejectionRates(nBlocks,jChain)>op.RejectionRateMax);
        nBelow = nBelow + (RejectionRates(nBlocks,jChain)<op.RejectionRateMin);
    end
    fprintf('\nNumber of chains with rejection rate too high: %.0f',nAbove)
    fprintf('\nNumber of chains with rejection rate too low: %.0f\n\n',nBelow)
    if nAbove==nBelow
        jConfirm = jConfirm+1;
        if jConfirm>op.NConfirm
            fprintf('Results confirmed!\n')
        else
            fprintf('Rejection rate in range. Confirming results...\n')
        end
    else
        lastChange = currChange;
        currChange = (nAbove<nBelow)-(nAbove>nBelow);
        fprintf('Direction: %.0f\n',currChange)
        if lastChange+currChange==0
            nRevertDirection = nRevertDirection+1;
            fprintf('Direction reverted %.0f times\n',nRevertDirection)
        end
        if nRevertDirection>op.MaxReverseDirection
            fprintf(['Direction reverted too many times. Changing ' ...
                     'increments...\n'])
            jIncrement = jIncrement+1;
            nRevertDirection = 0;
        end
        if jIncrement>length(op.ScaleIncrements)
            jConfirm = jConfirm + 1;
            if jConfirm>op.NConfirm
                fprintf('Results confirmed!\n')
            else
                fprintf(['Increments changed too many times. Confirming ',...
                         'results...\n'])
            end
        else
            jumpScaleOld = jumpScale;
            jumpScale = jumpScale + currChange*op.ScaleIncrements(jIncrement);
            if jumpScale<op.MinSearchScale
                jumpScale = op.MinSearchScale;
                fprintf(['Minimum jump scale breached. Setting it to minimum ' ...
                         'level: %.2f\n'],jumpScale)
                jConfirm = jConfirm + 1;
            else
                fprintf('Scale changed to %.2f\n',jumpScale)
                jConfirm = (jumpScaleOld==jumpScale);
            end
        end            
    end
end

%% Show rejection rates: 
% no longer showing this because it simply duplicates information already there
% fprintf('\nResults for Jump Scale calibration:\n\n')
% for jBlock=1:nBlocks
%     for jChain=1:op.NChains
%         fprintf(['Block %03.0f, Chain %02.0f: JumpScale = %4.2f, ',...
%                  'rejection rate = %4.1f%%\n'],...
%                 jBlock,jChain,JumpScaleFactors(jBlock),...
%                 RejectionRates(jBlock,jChain)*100)
%     end
%     fprintf('\n')
% end

%% show reason to stop
fprintf('\nReason to stop: ')
if jIncrement>length(op.ScaleIncrements)
    fprintf('Increments changed too many times.\n\n')
else
    fprintf('Results confirmed.\n\n')
end

%% save workspace
save(sprintf('%s-mcmc-%.0f-calibration',obj.Name,obj.Post.MCMCStage))

%% Clean up
if ~op.KeepFilesCalibrate
    for jChain=1:op.NChains
        delete([fn{jChain},'.log']);
        delete([fn{jChain},'.mat']);
    end
end
 



