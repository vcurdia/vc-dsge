function draw = dksmoother(mats,data,isDrawStates)

% dksmoother
% 
% Draw from the Durbin Koopman Disturbance Smoother.
%
% Usage:
%   draw = dksmoother(mats,data)
%   draw = dksmoother(mats,data,isDrawStates)
% 
% NOTE: assume that states have zero steady state
% 
% Options:
%   - isDrawStates (logical)
%     if set to 1, states are recursively drawn from conditional distributions.
%     if set to 0, states are set to their conditional distribution mean.
%     default: 1
%
% See also:
% DSGE.Model, DSGE.Model.states
%
% .............................................................................
% 
% Created: April 14, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Preamble

if nargin<3, isDrawStates = 1; end

sig00 = mats.KF.sig00;
s00 = mats.KF.s00;
G1 = mats.REE.G1;
G2 = mats.REE.G2;
H = mats.ObsEq.H;
idxNaN = isnan(data);

[nStateVar,nShockVar] = size(G2);
[T,nObsVar] = size(data);


%% Generate random states
Er = isDrawStates*normrnd(0,1,nShockVar,T);
Sr = zeros(nStateVar,T+1);
try 
    Sr(:,1) = s00+isDrawStates*mvnrnd(zeros(1,nStateVar),sig00)';
catch
%   fprintf('Warning: sig00 no semidefinite positive. using only diag elements.\n');
    sig00tmp = diag(sig00);
    sig00tmp(sig00tmp<0)=0;
    sig00 = diag(sig00tmp);
    Sr(:,1) = s00+isDrawStates*mvnrnd(zeros(1,nStateVar),sig00)';
end
for t=1:T
    Sr(:,t+1) = G1*Sr(:,t)+G2*Er(:,t);
end

%% Generate difference between observables and random series
DataDiff = data'-mats.KF.ObsVarBar-H*Sr(:,2:T+1);

%% Run KF on Xs
stt = s00;
sigtt = sig00;
r = zeros(nStateVar,T);
K = cell(T,1);
Om = G2*G2';
for t=1:T
    Ht = H(~idxNaN(t,:),:);
    sigtt1 = G1*sigtt*G1'+Om;
    Ft = Ht*sigtt1*Ht';
    vt = DataDiff(~idxNaN(t,:),t)-Ht*G1*stt;
    Kt = sigtt1*Ht'/Ft;
    stt = G1*stt+Kt*vt;
    sigtt = (eye(nStateVar)-Kt*Ht)*sigtt1;
    r(:,t) = Ht'*(Ft\vt);
    K{t} = Kt;
end

%% Run DK
for t=T-1:-1:1
    r(:,t) = r(:,t)+(eye(nStateVar)-...
                     H(~idxNaN(t,:),:)'*K{t}')*G1'*r(:,t+1);
end
r0 = G1'*r(:,1);

%% Get shocks and states
E = G2'*r;
S = zeros(nStateVar,T+1);
S(:,1) = s00+sig00*r0;
for t=1:T
  S(:,t+1) = G1*S(:,t)+G2*E(:,t);
end
StateDraw = Sr+S;
draw.StateVar = StateDraw(:,2:T+1);
draw.ShockVar = Er+E;
draw.StateVar0 = StateDraw(:,1);

