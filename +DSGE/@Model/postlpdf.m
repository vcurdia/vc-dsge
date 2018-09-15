function p = postlpdf(obj,x,varargin)

% postlpdf
% 
% posterior log-pdf function
%
% see also:
% DSGE.Posterior
%
% ............................................................................
%
% Created: March 23, 2017
% Copyright (C) 2017 Vasco Curdia

%% Options
op.FID = 1;
op.Verbose = 0;
op = updateoptions(op,varargin{:});

%% convert param vector
if size(x,1)<obj.Param.N
    x = obj.expandparam(x);
end

%% Evaluate prior log-pdf
p = obj.priorlpdf(x);
if p==-inf, return, end

%% Get Mats
Mats = obj.solveree(x,op);
if ~Mats.Status
    p = -inf;
    return
end

%% Kalman Filter
stt = Mats.KF.s00;
sigtt = Mats.KF.sig00;
StateVartt = zeros(obj.StateVar.N,obj.Data.T);
SIGtt = zeros(obj.StateVar.N,obj.StateVar.N,obj.Data.T);
DataDetrended = obj.Data.Values-Mats.KF.ObsVarBar';
for t=1:obj.Data.T
    idxNoNaN = ~isnan(DataDetrended(t,:));
    [stt,sigtt,lh,ObsVarhat] = kf(DataDetrended(t,idxNoNaN)',...
                                  Mats.ObsEq.H(idxNoNaN,:),stt,sigtt,...
                                  Mats.REE.G1,Mats.REE.G2);
    if t>obj.Data.NPreSample
        p = p + lh*[1;1];
    end
    StateVartt(:,t) = stt;
    SIGtt(:,:,t) = sigtt;
end

%% Add normalization
p = p - sum(sum(~isnan(obj.Data.Values(obj.Data.NPreSample+1:obj.Data.T,:))))...
    /2*log(2*pi);


