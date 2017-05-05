function op = Table(varargin)

% DSGE.Options.Table
%
% Set Table options
%
% Created: March 21, 2017 
% Copyright 2017 by Vasco Curdia

%% Default Options
op.Precision = 3;
op.MaxRows = 35;
op.MoveLeft = 1; 
op.Lines = [];

%% update options
op = updateoptions(op,varargin{:});
