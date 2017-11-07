function Mats = mats(obj,x,varargin)

% mats
% 
% Generates mats for a given parameter vector
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
op.NumsolveGuess = [];
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
    ObsEq.HBar = obj.MatFcn.Obseq.HBar(xAll);
    ObsEq.H = obj.MatFcn.Obseq.HB(xAll);
end
if op.StoreObsEq
    Mats.ObsEq = ObsEq;
end

HERE HERE HERE

%% State equation matrices
StateEq.GammaBar = [...
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    ];

StateEq.Gamma0 = [...
    -1, 0, 0, -1/phigamma, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, beta*etagamma, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 1/omega, 1/omega, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -beta*etagamma*phigamma, (beta*etagamma)/(beta*etagamma - 1), -beta*etagamma*phigamma, 0, 0, 0;
    0, 0, -beta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    ];

StateEq.Gamma1 = [...
    -1, 0, 0, 0, -1/phigamma, 0, 1/phigamma, 0, 0, 0, 0, 0, 0;
    -1, 0, 0, 0, 0, beta*etagamma^2 + 1, 0, 0, 0, 0, 0, -etagamma, etagamma;
    0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1/omega, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, - omega - phigamma*(beta*etagamma^2 + 1), 0, -etagamma*phigamma, 0, 0, etagamma*phigamma;
    phigamma*xi, 0, -1, 0, 0, omega*xi, 0, 0, 0, 0, 1, 0, 0;
    0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -phipi*(rho - 1), -1, -(phix*(rho - 1))/4, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    ];

StateEq.Gamma4 = [...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -zeta, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, rho, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, rhodelta, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, rhogamma, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rhou, 0, 0;
    0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    ];

StateEq.Gamma2 = [...
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, sigmai/400;
    sigmadelta/400, 0, 0, 0;
    0, sigmagamma/400, 0, 0;
    0, 0, sigmau/400, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    ];

StateEq.Gamma3 = eye(13);

cv = (all(StateEq.Gamma0(1:13,:)==0,2)~=0);
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

% Solve REE
if op.SolveREE
    [REE,fmat,fwt,ywt,gev] = solveree(StateEq,...
        op.GensysAuthor,op.fid,op.verbose,op.gensys{:});
    Mats.REE = REE;
    if ~all(REE.eu==1);
        Mats.Status = 0;
        Mats.StatusMessage = [Mats.StatusMessage,'REE solution not normal.'];
    end
end

% Kalman Filter matrices
if op.SolveREE && op.StoreKF
    if all(Mats.REE.GBar(:)==0)
        KF.StateVarBar = zeros(13,1);
    else
        KF.StateVarBar = (eye(13)-REE.G1)\REE.GBar;
    end
    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;

    KF.s00 = zeros(13,1);

    [sig00,sig00rc] = lyapcsd(REE.G1,REE.G2*REE.G2');
    sig00 = real(sig00); sig00 = (sig00+sig00')/2;
    if sig00rc~=0
        Mats.Status = 0;
        Mats.StatusMessage = [Mats.StatusMessage,'Could not find unconditional variance.'];
        if op.verbose
            fprintf(fid,'Warning: Could not find unconditional variance.\n');
        end

    end

    KF.sig00 = sig00;
    KF.sig00rc = sig00rc;
    Mats.KF = KF;
end

% Auxiliary equations matrices
if op.StoreAuxEq || op.StoreAuxREE
    AuxEq.PhiBar = [...
        0;
        ];

    AuxEq.Phi1 = [...
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
        ];

    AuxEq.Phi2 = [...
        0, 0, 0, 0;
        ];

    AuxEq.Phi3 = [...
        0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        ];

    AuxEq.Phi4 = [...
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        ];

    if op.StoreAuxEq
        Mats.AuxEq = AuxEq;
    end
    if op.SolveREE && op.StoreAuxREE
        if ~isempty(REE.G1)
            AuxREE.GBar = AuxEq.PhiBar+AuxEq.Phi3*REE.GBar+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.GBar;
            AuxREE.G1 = AuxEq.Phi4+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G1;
            AuxREE.G2 = AuxEq.Phi2+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G2;
        else
            AuxREE.GBar = [];
            AuxREE.G1 = [];
            AuxREE.G2 = [];
        end
        Mats.AuxREE = AuxREE;
    end
end
end
