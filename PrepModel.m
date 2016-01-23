function s = PrepModel(s,Model)

% PrepModel
%
% Analyzes the model and generates code to evaluate the model for a given
% parameter vector.
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
%             x_ss refers to steady state of x(t)
%
% See also:
% RunMyDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

Action = 'PrepModel';

% check if model already prepared
if isfield(s.Status,Action) && s.Status.(Action), return, end

fprintf('Analyzing DSGE model...\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% -------------------------------------------------------------------


%% -------------------------------------------------------------------

%% Finish up

% Record Model as prepared
s.Status.(Action) = 1;

% Time Elapsed
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
