function Mats = evalmats(obj,x,varargin)

% evalmats
% 
% Evaluates DSGE mats for a given parameter vector
% 
% See also:
% DSGE, SetupMyDSGE, Model
%
% .............................................................................
% 
% Created: November 6, 2017 by Vasco Curdia
% 
% Copyright 2016-2017 by Vasco Curdia


%% Options
op.StoreStateEq = 1;
op.StoreObsEq = 1;
op.StoreKF = 1;
op.StoreAuxEq = 1;
op.StoreAuxREE = 1;
op.SolveREE = 1;
op.fid = 1;
op.verbose = 0;
op.GensysAuthor = 'CS';
op.gensys = {};
op.NumSolveGuess = [];
op = updateoptions(op,varargin{:});

%% Initiate Status
Mats.Status = 1;
Mats.StatusMessage = '';
Mats.Param = x;

%% Auxiliary parameters
if isempty(op.NumSolveGuess)
    xns0 = ones(obj.NumSolveParam.N,1);
end
NumSolveOptions = optimoptions(@fsolve);
NumSolveOptions.Display = 'off';
[xns1,NumSolveResidual,NumSolveRC,NumSolveOutput] = ...
    fsolve(@(xns)obj.MatFcn.NumSolveEq([x;xns]),xns0,NumSolveOptions);
Mats.AuxParam = [xns1;obj.MatFcn.CompoundParam([x;xns1])];
Mats.NumSolveParamRC = NumSolveRC;
if NumSolveRC~=1
    Mats.Status = 0;
    Mats.StatusMessage = [Mats.StatusMessage,...
                        'NumSolveParam solution not normal.'];
    if op.verbose
        fprintf(fid,'Warning: NumSolveParam solution not normal.\n');
    end
end
xAll = [Mats.Param;Mats.AuxParam];

%% Observation equations
if op.StoreObsEq || op.StoreKF
    ObsEq.HBar = obj.MatFcn.ObsEq.HBar(xAll);
    ObsEq.H = obj.MatFcn.ObsEq.H(xAll);
end
if op.StoreObsEq
    Mats.ObsEq = ObsEq;
end

%% State equations
StateEq.GammaBar = obj.MatFcn.StateEq.GammaBar(xAll);
StateEq.Gamma0 = obj.MatFcn.StateEq.Gamma0(xAll);
StateEq.Gamma1 = obj.MatFcn.StateEq.Gamma1(xAll);
StateEq.Gamma2 = obj.MatFcn.StateEq.Gamma2(xAll);
StateEq.Gamma3 = eye(obj.StateVar.N);
StateEq.Gamma4 = obj.MatFcn.StateEq.Gamma4(xAll);
cv = (all(StateEq.Gamma0(1:obj.StateVar.N,:)==0,2)~=0);
StateEq.Gamma0(cv,:) = -StateEq.Gamma1(cv,:);
StateEq.Gamma1(cv,:) = StateEq.Gamma4(cv,:);
StateEq.Gamma3(:,cv) = [];
if ~all(all(StateEq.Gamma4(~cv,:)==0,2))
    error('Incorrect system reduction')
end
StateEq = rmfield(StateEq,'Gamma4');
if op.StoreStateEq
    Mats.StateEq = StateEq;
end

%% Solve REE
if op.SolveREE
    [REE,fmat,fwt,ywt,gev] = solveree(StateEq,...
        op.GensysAuthor,op.fid,op.verbose,op.gensys{:});
    Mats.REE = REE;
    if ~all(REE.eu==1);
        Mats.Status = 0;
        Mats.StatusMessage = [Mats.StatusMessage,'REE solution not normal.'];
    end
end

%% Kalman Filter matrices
if op.SolveREE && op.StoreKF
    if all(Mats.REE.GBar(:)==0)
        KF.StateVarBar = zeros(obj.StateVar.N,1);
    else
        KF.StateVarBar = (eye(obj.StateVar.N)-REE.G1)\REE.GBar;
    end
    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;
    if isempty(obj.KFInitState)
        KF.s00 = zeros(obj.StateVar.N,1);
    else
        KF.s00 = obj.KFInitState;
    end
    if isempty(obj.KFInitVariance)
        [sig00,sig00rc] = lyapcsd(REE.G1,REE.G2*REE.G2');
        sig00 = real(sig00); 
        sig00 = (sig00+sig00')/2;
        if sig00rc~=0
            Mats.Status = 0;
            Mats.StatusMessage = [Mats.StatusMessage, ...
                                'Could not find unconditional variance.'];
            if op.verbose
                fprintf(fid, ...
                        'Warning: Could not find unconditional variance.\n');
            end
        end
    else
        sig00 = obj.KFInitVariance;
        sig00rc = 0;
    end
    KF.sig00 = sig00;
    KF.sig00rc = sig00rc;
    Mats.KF = KF;
end

%% Auxiliary equations matrices
if op.StoreAuxEq || op.StoreAuxREE
    AuxEq.PhiBar = obj.MatFcn.AuxEq.PhiBar(xAll);
    AuxEq.Phi1 = obj.MatFcn.AuxEq.Phi1(xAll);
    AuxEq.Phi2 = obj.MatFcn.AuxEq.Phi2(xAll);
    AuxEq.Phi3 = obj.MatFcn.AuxEq.Phi3(xAll);
    AuxEq.Phi4 = obj.MatFcn.AuxEq.Phi4(xAll);
    if op.StoreAuxEq
        Mats.AuxEq = AuxEq;
    end
    if op.SolveREE && op.StoreAuxREE
        if ~isempty(REE.G1)
            AuxREE.GBar = AuxEq.PhiBar + AuxEq.Phi3*REE.GBar ...
                + (AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.GBar;
            AuxREE.G1 = AuxEq.Phi4 + (AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G1;
            AuxREE.G2 = AuxEq.Phi2 + (AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G2;
        else
            AuxREE.GBar = [];
            AuxREE.G1 = [];
            AuxREE.G2 = [];
        end
        Mats.AuxREE = AuxREE;
    end
end
