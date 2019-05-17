function Mats = solveree(obj,x,varargin)

% solveree
%
% Evaluates model matrices and useses gensys to solve for the REE and KF
%
% Option: FastGensys (logical)
% If set to 0 (default) it uses the original Chris Sims verion. If set to 1 it
% uses the fast gensys from Jae Won and if it does not yield a normal solution 
% it runs the original gensys.
% 
% Other options available
%
% Usage:
%   Mats = solveree(obj,x)
%   Mats = solveree(obj,x,varargin)
%
% See also:
% DSGE.Model, DSGE.Model.genmats, gensys, gensysvb, fastgensysJaeWonvb
%
% ...........................................................................
% 
% Created: January 25, 2016 by Vasco Curdia
% Copyright 2016-2018 by Vasco Curdia


% Default options
op.FID = 1;
op.Verbose = 0;
op.FastGensys = 0;
op.Div = [];
op.RealSmall = [];
op.UsePinv = 0;
op = updateoptions(op,obj.GensysOptions,varargin{:});

%% set fid
fid = op.FID;

%% evaluate model mats
Mats = obj.mats(x,op);
if ~Mats.Status
    if op.Verbose
        fprintf(fid,'Warning: %s\n',Mats.StatusMessage);
    end
    return
end
    
gsInput = [Mats.StateEq.Gamma0(:);Mats.StateEq.Gamma1(:);...
           Mats.StateEq.GammaBar(:);Mats.StateEq.Gamma2(:);...
                           Mats.StateEq.Gamma3(:)];
if any(isnan(gsInput)) || any(isinf(gsInput))
    Mats.Status = 0;
    Mats.StatusMessage = 'Model equations contain NAN or INF elements.';
    if op.Verbose
        fprintf(fid,'Warning: %s\n',Mats.StatusMessage);
    end
    return
end


%% Run Gensys
REE.GBar = [];
REE.G1 = [];
REE.G2 = [];
REE.eu = [0;0];
if op.FastGensys
    [REE.G1,REE.GBar,REE.G2,fmat,fwt,ywt,gev,REE.eu] = ...
        fastgensysJaeWonvb(Mats.StateEq.Gamma0,Mats.StateEq.Gamma1,...
                           Mats.StateEq.GammaBar,Mats.StateEq.Gamma2,...
                           Mats.StateEq.Gamma3,...
                           op.FID,op.Verbose,op.Div,op.RealSmall,op.UsePinv);
end
if ~op.FastGensys || ~all(REE.eu(:)==1)
    [REE.G1,REE.GBar,REE.G2,fmat,fwt,ywt,gev,REE.eu] = ...
        gensysvb(Mats.StateEq.Gamma0,Mats.StateEq.Gamma1,...
                 Mats.StateEq.GammaBar,Mats.StateEq.Gamma2,...
                 Mats.StateEq.Gamma3,...
                 op.FID,op.Verbose,op.Div,op.RealSmall,op.UsePinv);
end
Mats.REE = REE;
if ~all(REE.eu==1);
    Mats.Status = 0;
    Mats.StatusMessage = 'REE solution not normal.';
    if op.Verbose
        fprintf(fid,'Warning: %s\n',Mats.StatusMessage);
    end
    return
end

%% Kalman Filter matrices
if obj.ObsVar.N>0
    if all(REE.GBar(:)==0)
        KF.StateVarBar = zeros(obj.StateVar.N,1);
    else
        KF.StateVarBar = (eye(obj.StateVar.N)-REE.G1)\REE.GBar;
    end
    KF.ObsVarBar = Mats.ObsEq.HBar + Mats.ObsEq.H*KF.StateVarBar;

    if ~isempty(obj.KFInitState)
        s00 = obj.KFInitState;
    else
        KF.s00 = zeros(obj.StateVar.N,1);
    end

    if ~isempty(obj.KFInitVariance)
        sig00 = obj.KFInitVariance;
        sig00rc = 0;
    else
        [sig00,sig00rc] = lyapcsdvb(REE.G1,REE.G2*REE.G2',op.Verbose);
        sig00 = real(sig00); 
        sig00 = (sig00+sig00')/2;
        if sig00rc~=0
            Mats.Status = 0;
            Mats.StatusMessage = 'Could not find unconditional variance.';
            if op.Verbose
                fprintf(fid,'Warning: %s\n',Mats.StatusMessage);
            end
        end
    end
    KF.sig00 = sig00;
    KF.sig00rc = sig00rc;
    Mats.KF = KF;
end
