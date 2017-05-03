function sim(obj,varargin)

% sim
% 
% Simulate DSGE model
% 
% See also:
% DSGE.Model, setupMyDSGE, runDSGE
%
% .............................................................................
% 
% Created: May 1, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Default Options
op.FNSuffix = '';
op.Prior = [];
op.Post = [];
op.Data = [];

op = updateoptions(op,varargin{:});

%% run sim
obj.irf(xd,op)
obj.vd(xd,op)
obj.states(data,xd,op)
obj.sd(data,xd,op)

