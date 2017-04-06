function op = Figure(varargin)

% DSGE.Options.Figure
%
% Set Figure options
%
% ...........................................................................
% 
% Created: March 21, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%% Default Options
op.Visible = 'off';

%% update options
op = updateoptions(op,varargin{:});
