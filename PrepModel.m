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

%% Prepare variables

% Check for non-specified model fields
if ~isfield(Model,'AuxParam'), Model.AuxParam = cell(0,2); end
if ~isfield(Model,'AuxVar'), Model.AuxVar = cell(0,2); end
s.Model = Model;

% Parameters
Param.Names = {Model.Param{:,1}}';
if size(Param,2)==4
    Param.PrettyNames = Param.Names;
else
    Param.PrettyNames = {Model.Param{:,5}}';
end
Param.PriorDist = {Model.Param{:,2}}';
Param.PriorMean = [Model.Param{:,3}]';
Param.PriorSE = [Model.Param{:,4}]';
s.Param = Param;
s.nParam = length(Param.Names);

% Auxiliary Parameters
AuxParam.Names = {Model.AuxParam{:,1}}';
if size(AuxParam,2)==2
    AuxParam.PrettyNames = AuxParam.Names;
else
    AuxParam.PrettyNames = {Model.AuxParam{:,3}}';
end
s.AuxParam = AuxParam;
s.nAuxParam = length(AuxParam.Names);

% Observation variables
s.nObsVar = length(Model.ObsVar);
ObsVar_t = sym(zeros(1,s.nObsVar));
for j=1:s.nObsVar
    eval(sprintf('syms %s_t',Model.ObsVar{j}))
    ObsVar_t(j) = eval([Model.ObsVar{j},'_t']);
end

% State space variables
s.nStateVar = length(Model.StateVar);
StateVar_t = sym(zeros(1,s.nStateVar)); 
StateVar_tF = sym(zeros(1,s.nStateVar)); 
StateVar_tL = sym(zeros(1,s.nStateVar));
for j=1:s.nStateVar
    eval(sprintf('syms %1$s_t %1$s_tF %1$s_tL',Model.StateVar{j}))
    StateVar_t(j) = eval([Model.StateVar{j},'_t']);
    StateVar_tF(j) = eval([Model.StateVar{j},'_tF']);
    StateVar_tL(j) = eval([Model.StateVar{j},'_tL']);
end

% Shocks
s.nShockVar = length(Model.ShockVar);
ShockVar_t = sym(zeros(1,s.nShockVar));
for j=1:s.nShockVar
    eval(['syms ',Model.ShockVar{j},'_t'])
    ShockVar_t(j) = eval([Model.ShockVar{j},'_t']);
end

% constant variable
syms one

% AuxVar
s.nAuxVar = size(Model.AuxVar,1);
for j=1:s.nAuxVar
    eval([Model.AuxVar{j,1},'_t = sym(',Model.AuxVar{j,2},');']);
end
check for tF var in definition. If not, create tF variable.
check for tL var in definition. If not, create tL variable.
    
keyboard

%% -------------------------------------------------------------------

%% Finish up

% Record Model as prepared
s.Status.(Action) = 1;

% Time Elapsed
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
