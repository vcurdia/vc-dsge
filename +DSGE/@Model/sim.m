function sim(obj,xd,varargin)

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
op.List = {'IRF','VD','States','SD'};
op.Data = [];

op = updateoptions(op,varargin{:});

%% run sim
if ismember('IRF',op.List), obj.irf(xd,op); end
if ismember('VD',op.List), obj.vd(xd,op), end
if ~isempty(op.Data)
    if ismember('States',op.List), obj.states(xd,op), end
    if ismember('SD',op.List), obj.sd(xd,op), end
end

