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
% DSGESetup
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

% check if model already prepared
if s.Done.PrepModel, return, end

fprintf('Analyzing DSGE model...\n')

% Set Timer
s.TimeElapsed.PrepModel = toc();

%% -------------------------------------------------------------------


%% -------------------------------------------------------------------

%% Finish up

% Record Model as prepared
s.Done.PrepModel = 1;

% Time Elapsed
s.TimeElapsed.PrepModel = toc-TimeElapsed.GenSymVars;

%% -------------------------------------------------------------------
