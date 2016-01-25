function s = PrepModel(s,Model,GensysAuthor)

% PrepModel
%
% Analyzes the model and generates code to evaluate the model for a given
% parameter vector.
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
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

fprintf('\n*** Analyzing DSGE model\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% -------------------------------------------------------------------

%% Prepare variables and equations

fprintf('Generating symbolic variables and systems of equations...\n')

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
for j=1:s.nParam
    eval(['syms ',Param.Names{j}]);
end

% Auxiliary Parameters
AuxParam.Names = {Model.AuxParam{:,1}}';
if size(AuxParam,2)==2
    AuxParam.PrettyNames = AuxParam.Names;
else
    AuxParam.PrettyNames = {Model.AuxParam{:,3}}';
end
s.AuxParam = AuxParam;
s.nAuxParam = length(AuxParam.Names);
for j=1:s.nAuxParam
    eval(['syms ',AuxParam.Names{j}]);
end

% Observation variables
s.nObsVar = length(Model.ObsVar);
ObsVar_t = sym(zeros(1,s.nObsVar));
for j=1:s.nObsVar
    jv = [Model.ObsVar{j},'_t'];
    eval(['syms ',jv]);
    ObsVar_t(j) = eval(jv);
end

% State space variables
s.nStateVar = length(Model.StateVar);
StateVar_t = sym(zeros(1,s.nStateVar)); 
StateVar_tF = sym(zeros(1,s.nStateVar)); 
StateVar_tL = sym(zeros(1,s.nStateVar));
for j=1:s.nStateVar
    jv = [Model.StateVar{j},'_t'];
    eval(sprintf('syms %1$s %1$sF %1$sL',jv))
    StateVar_t(j) = eval(jv);
    StateVar_tF(j) = eval([jv,'F']);
    StateVar_tL(j) = eval([jv,'L']);
end

% Shocks
s.nShockVar = length(Model.ShockVar);
ShockVar_t = sym(zeros(1,s.nShockVar));
for j=1:s.nShockVar
    jv = [Model.ShockVar{j},'_t'];
    eval(['syms ',jv]);
    ShockVar_t(j) = eval(jv);
end

% constant variable
syms one

% Auxiliary variables
s.nAuxVar = size(Model.AuxVar,1);
AuxVar_t = sym(zeros(s.nAuxVar,1));
for j=1:s.nAuxVar
    jv = [Model.AuxVar{j,1},'_t'];
    eval([jv,' = ',Model.AuxVar{j,2},';'])
    AuxVar_t(j) = eval(jv);
    % if the expression has no leads then can define a lead for it
    if all(jacobian(AuxVar_t(j),StateVar_tF)==0)
        eval([jv,'F = subs(',jv,',[StateVar_t,StateVar_tL]',...
              ',[StateVar_tF,StateVar_t]);'])
        % if expreassion has no leads or lags then can define lag for it. 
        % notice that it does not make sense to define a lag if there are 
        % leads in it, and that's why the check for lags is inside the check 
        % for leads
        if all(jacobian(AuxVar_t(j),StateVar_tL)==0)
            eval([jv,'L = subs(',jv,',[StateVar_tF,StateVar_t]',...
                  ',[StateVar_t,StateVar_tL]);'])
        end
    end
end

% Build Observation equations
s.nObsEq = length(Model.ObsEq);
if s.nObsEq~=s.nObsVar
    error(['Number of observables (%\.0f) is different from number of ' ...
           'observation equations (%.0f).'],s.nObsVar,s.nObsEq)
end
ObsEq = sym(zeros(s.nObsEq,1));
for j=1:s.nObsEq
    ObsEq(j) = sym(Model.ObsEq{j});
end

% Build State equations
s.nStateEq = length(Model.StateEq);
if s.nStateEq~=s.nStateVar
    error(['Number of state variables (%\.0f) is different from number of ' ...
           'state equations (%.0f).'],s.nStateVar,s.nStateEq)
end
StateEq = sym(zeros(s.nStateEq,1));
for j=1:s.nStateEq
    StateEq(j) = sym(Model.StateEq{j});
end


%% Generate matrices and MakeMats

if ~exist('GensysAuthor','var'),GensysAuthor='CS';end

fprintf('Generating Mats for model evaluation\n')
s.FileName.Mats = sprintf('%sMats',s.Spec);

% Initiate file
fidMats = fopen([s.FileName.Mats,'.m'],'wt');
fprintf(fidMats,...
        'function Mats=%s(x,fid,verbose,varargin)\n\n',...
        s.FileName.Mats);
fprintf(fidMats,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs \n\n',clock);

% close file
fclose(fidMats);

    


%% -------------------------------------------------------------------

%% Finish up

% Record Model as prepared
s.Status.(Action) = 1;

% Time Elapsed
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
