function sim(obj,data,xd,varargin)

% sim
% 
% Simulate DSGE model
% 
% See also:
% DSGE, SetupMyDSGE
%
% .............................................................................
% 
% Created: May 1, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Default Options
op.FNSuffix = '';

op = updateoptions(op,varargin{:});

%% run sim
obj.irf(xd,op)
obj.vd(xd,op)
obj.states(data,xd,op)
obj.sd(data,xd,op)

