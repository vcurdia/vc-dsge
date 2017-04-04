function jumpScale = calibratejump(obj,varargin)

% calibratejump
% 
% Calibrate MCMC jump distribution
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.SampleID = length(obj.MCMCSample)+1;
op.NChains = 4;
op.JumpScale = 2.4;
op.NDrawsCalibrateJump = 1000;
op.KeepFilesCalibrateJump = 0; 
op.ScaleIncrements = [0.2,0.05,0.01];
op.NConfirm = 1;
op.RejectionRateMax = 0.85;
op.RejectionRateMin = 0.70; 
op.MaxReverseDirection = 2; 
op.MinSearchScale = 0.1; 

op = updateoptions(op,varargin{:});

%% Preparations

sid = op.SampleID;
pIdx = obj.EstimateIdx;

fprintf('\n*** Calibrating Jump Distribution for MCMC Sample %.0f\n',sid)
ttName = sprintf('CalibrateJumpMCMCSample%.0f',sid);
obj.TimeElapsed.start(ttName)

jumpScale = op.JumpScale;

op.Chain.Augment = 0;
op.Chain.NDraws = op.NDrawsCalibrateJump;
op.Chain.x0 = [];
op.Chain.NRejections = 0;

fn = cell(op.NChains,1);
for jChain=1:op.NChains
    fn{jChain} = sprintf('%s_MCMC_Sample_%.0f_Chain_%.0f_CalibrateJump',...
                         obj.Model.Name,sid,jChain);
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
    op.Chain.JumpVar = jumpScale^2/obj.NEstimate*obj.Var(pIdx,pIdx);
    RejectionRates(nBlocks,:) = zeros(1,op.NChains);
    parfor jChain=1:op.NChains
        opj = op.Chain;
        opj.fn = fn{jChain};
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
            jConfirm = jConfirm+1;
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
            else
                fprintf('Scale changed to %.2f\n',jumpScale)
            end
            jConfirm = (jumpScaleOld==jumpScale);
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


%% Clean up
if ~op.KeepFilesCalibrateJump
    for jChain=1:op.NChains
        delete([fn{jChain},'.log']);
        delete([fn{jChain},'.mat']);
    end
end
 

%% Finish up
obj.TimeElapsed.stop(ttName)


