function s = PrepModel(s)

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
% Created: January 27, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

Action = 'PrepModel';

% % check if model already prepared
% if isfield(s.Status,Action) && s.Status.(Action), return, end

fprintf('\n*** Analyzing DSGE model\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

% Check settings
if ~isfield(s,'GensysAuthor'), s.GensysAuthor = 'CS';end


%% -------------------------------------------------------------------

%% Prepare variables and equations

fprintf('Generating symbolic variables and systems of equations...\n')

% Check for non-specified model fields
if ~isfield(s,'AuxParam'), s.AuxParam = cell(0,2); end
if ~isfield(s,'AuxVar'), s.AuxVar = cell(0,2); end

% Parameters
Param.Names = {s.Param{:,1}}';
if size(s.Param,2)==4
    Param.PrettyNames = Param.Names;
else
    Param.PrettyNames = {s.Param{:,5}}';
end
Param.PriorDist = {s.Param{:,2}}';
Param.PriorMean = [s.Param{:,3}]';
Param.PriorSE = [s.Param{:,4}]';
s.Param = Param;
s.nParam = length(Param.Names);
for j=1:s.nParam
    eval(['syms ',Param.Names{j}]);
end

% Auxiliary Parameters
AuxParam.Names = {s.AuxParam{:,1}}';
if size(s.AuxParam,2)==2
    AuxParam.PrettyNames = AuxParam.Names;
else
    AuxParam.PrettyNames = {s.AuxParam{:,3}}';
end
AuxParam.Expressions = {s.AuxParam{:,2}}';
s.AuxParam = AuxParam;
s.nAuxParam = length(AuxParam.Names);
for j=1:s.nAuxParam
    eval(['syms ',AuxParam.Names{j}]);
end

% Observation variables
s.nObsVar = length(s.ObsVar);
ObsVar_t = sym(zeros(1,s.nObsVar));
for j=1:s.nObsVar
    jv = [s.ObsVar{j},'_t'];
    eval(['syms ',jv]);
    ObsVar_t(j) = eval(jv);
end

% State space variables
s.nStateVar = length(s.StateVar);
StateVar_t = sym(zeros(1,s.nStateVar)); 
StateVar_tF = sym(zeros(1,s.nStateVar)); 
StateVar_tL = sym(zeros(1,s.nStateVar));
for j=1:s.nStateVar
    jv = [s.StateVar{j},'_t'];
    eval(sprintf('syms %1$s %1$sF %1$sL',jv))
    StateVar_t(j) = eval(jv);
    StateVar_tF(j) = eval([jv,'F']);
    StateVar_tL(j) = eval([jv,'L']);
end

% Shocks
s.nShockVar = length(s.ShockVar);
ShockVar_t = sym(zeros(1,s.nShockVar));
for j=1:s.nShockVar
    jv = [s.ShockVar{j},'_t'];
    eval(['syms ',jv]);
    ShockVar_t(j) = eval(jv);
end

% constant variable
syms one

% Auxiliary variables
s.AuxEq = {s.AuxVar{:,2}}.';
s.AuxVar = {s.AuxVar{:,1}}.';
s.nAuxVar = length(s.AuxVar);
AuxEq = sym(zeros(s.nAuxVar,1));
for j=1:s.nAuxVar
    jv = [s.AuxVar{j},'_t'];
    eval([jv,' = ',s.AuxEq{j},';'])
    AuxEq(j) = eval(jv);
    % if the expression has no leads then can define a lead for it
    if all(jacobian(AuxEq(j),StateVar_tF)==0)
        eval([jv,'F = subs(',jv,',[StateVar_t,StateVar_tL]',...
              ',[StateVar_tF,StateVar_t]);'])
        % if expreassion has no leads or lags then can define lag for it. 
        % notice that it does not make sense to define a lag if there are 
        % leads in it, and that's why the check for lags is inside the check 
        % for leads
        if all(jacobian(AuxEq(j),StateVar_tL)==0)
            eval([jv,'L = subs(',jv,',[StateVar_tF,StateVar_t]',...
                  ',[StateVar_t,StateVar_tL]);'])
        end
    end
end

% Build Observation equations
s.nObsEq = length(s.ObsEq);
if s.nObsEq~=s.nObsVar
    error(['Number of observables (%\.0f) is different from number of ' ...
           'observation equations (%.0f).'],s.nObsVar,s.nObsEq)
end
ObsEq = sym(zeros(s.nObsEq,1));
for j=1:s.nObsEq
    ObsEq(j) = sym(s.ObsEq{j});
end

% Build State equations
s.nStateEq = length(s.StateEq);
if s.nStateEq~=s.nStateVar
    error(['Number of state variables (%\.0f) is different from number of ' ...
           'state equations (%.0f).'],s.nStateVar,s.nStateEq)
end
StateEq = sym(zeros(s.nStateEq,1));
for j=1:s.nStateEq
    StateEq(j) = sym(s.StateEq{j});
end


%% Generate matrices

fprintf('Generating Mats for model evaluation\n')
s.FileName.Mats = sprintf('%sMats',s.Spec);

% Initiate file
fidMats = fopen([s.FileName.Mats,'.m'],'wt');
fprintf(fidMats,...
        'function Mats = %s(x,varargin)\n\n',...
        s.FileName.Mats);
fprintf(fidMats,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs\n',clock);


fprintf(fidMats,'\n%% Default options\n');
fprintf(fidMats,'op.StoreParam = 1;\n');
fprintf(fidMats,'op.StoreStateEq = 1;\n');
fprintf(fidMats,'op.StoreObsEq = 1;\n');
fprintf(fidMats,'op.StoreKF = 1;\n');
fprintf(fidMats,'op.StoreAuxEq = 1;\n');
fprintf(fidMats,'op.SolveREE = 1;\n');
fprintf(fidMats,'op.fid = 1;\n');
fprintf(fidMats,'op.verbose = 0;\n');

fprintf(fidMats,'\n%% Update options\n');
fprintf(fidMats,'if length(varargin)>0 && isstruct(varargin{1})\n');
fprintf(fidMats,'    op = varargin{1};\n');
fprintf(fidMats,'    varargin(1) = [];\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'for jop=1:(length(varargin)/2)\n'); 
fprintf(fidMats,'    op.(varargin{(jop-1)*2+1}) = varargin{jop*2};\n');
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Verify options\n');
fprintf(fidMats,'if op.StoreKF, op.SolveREE = 1; end\n');

fprintf(fidMats,'\n%% Map parameters\n');
for j=1:s.nParam
    fprintf(fidMats,'%s = x(%.0f);\n',Param.Names{j},j);
end
fprintf(fidMats,'if op.StoreParam\n');
for j=1:s.nParam
    fprintf(fidMats,'    Mats.Param.%s = x(%.0f);\n',Param.Names{j},j);
end
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Map auxiliary parameters\n');
for j=1:s.nAuxParam
    fprintf(fidMats,'%s = %s;\n',AuxParam.Names{j},AuxParam.Expressions{j});
end
fprintf(fidMats,'if op.StoreParam\n');
for j=1:s.nAuxParam
    fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',AuxParam.Names{j});
end
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Observation equations\n');
H0 = -jacobian(ObsEq,ObsVar_t);
SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
SymMats.ObsEq.H = H0\jacobian(ObsEq,StateVar_t);
MatNames = fieldnames(SymMats.ObsEq);
nCols = [1,s.nStateVar];
fprintf(fidMats,'if op.StoreObsEq || op.StoreKF\n');
for jM=1:length(MatNames)
    fprintf(fidMats,'    ObsEq.%s = [...\n',MatNames{jM});
    for jeq=1:s.nObsVar
        fprintf(fidMats,'       ');
        for jc=1:nCols(jM)
            fprintf(fidMats,' %s',...
                    char(eval(sprintf('SymMats.ObsEq.%s(jeq,jc)',...
                                      MatNames{jM}))));
            if jc==nCols(jM)
                fprintf(fidMats,';\n');
            else
                fprintf(fidMats,',');
            end
        end
    end
    fprintf(fidMats,'        ];\n');
end
fprintf(fidMats,'end\n');
fprintf(fidMats,'if op.StoreObsEq\n');
fprintf(fidMats,'    Mats.ObsEq = ObsEq;\n');
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% State equation matrices\n');
SymMats.StateEq.GammaBar = jacobian(StateEq,one);
SymMats.StateEq.Gamma0 = -jacobian(StateEq,StateVar_tF);
SymMats.StateEq.Gamma1 = jacobian(StateEq,StateVar_t);
SymMats.StateEq.Gamma4 = jacobian(StateEq,StateVar_tL);
SymMats.StateEq.Gamma2 = jacobian(StateEq,ShockVar_t);
MatNames = fieldnames(SymMats.StateEq);
nCols = [1,s.nStateVar,s.nStateVar,s.nStateVar,s.nShockVar];
for jM=1:length(MatNames)
    fprintf(fidMats,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:s.nStateVar
        fprintf(fidMats,'   ');
        for jc=1:nCols(jM)
            fprintf(fidMats,' %s',...
                    char(eval(sprintf('SymMats.StateEq.%s(jeq,jc)',...
                                      MatNames{jM}))));
            if jc==nCols(jM)
                fprintf(fidMats,';\n');
            else
                fprintf(fidMats,',');
            end
        end
    end
    fprintf(fidMats,'    ];\n\n');
end
fprintf(fidMats,'StateEq.Gamma3 = eye(%.0f);\n\n',s.nStateVar);
fprintf(fidMats,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',s.nStateVar);
fprintf(fidMats,'StateEq.Gamma0(cv,:) = -StateEq.Gamma1(cv,:);\n');
fprintf(fidMats,'StateEq.Gamma1(cv,:) = StateEq.Gamma4(cv,:);\n');
fprintf(fidMats,'StateEq.Gamma3(:,cv) = [];\n');
fprintf(fidMats,'if ~all(all(StateEq.Gamma4(~cv,:)==0,2))\n');
fprintf(fidMats,'    error(''Incorrect system reduction'')\n');
fprintf(fidMats,'end\n\n');
fprintf(fidMats,'StateEq = rmfield(StateEq,''Gamma4'');\n');
fprintf(fidMats,'if op.StoreStateEq\n');
fprintf(fidMats,'    Mats.StateEq = StateEq;\n');
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Solve REE\n');
fprintf(fidMats,'if op.SolveREE\n');
fprintf(fidMats,...
        '    [REE,fmat,fwt,ywt,gev] = SolveREE(StateEq,...\n');
fprintf(fidMats,...
        '        ''%s'',op.fid,op.verbose,varargin{:});\n',s.GensysAuthor);
fprintf(fidMats,'    Mats.REE = REE;\n');
fprintf(fidMats,'end\n');


fprintf(fidMats,'\n%% Kalman Filter matrices\n');
fprintf(fidMats,'if op.StoreKF\n');
fprintf(fidMats,'    if all(Mats.REE.GBar(:)==0)\n');
fprintf(fidMats,'        KF.StateVarBar = zeros(%.0f,1);\n',s.nStateVar);
fprintf(fidMats,'    else\n');
fprintf(fidMats,'        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
        s.nStateVar);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

if isfield(s,'KFinit') && isfield(s.KFinit,'State')
    fprintf(fidMats,'    s00 = [...\n');
    for jeq=1:s.nStateVar
        fprintf(fidMats,'        %.16f;\n',s.KFinit.State(jeq));
    end
    fprintf(fidMats,'        ];\n\n');
else
    fprintf(fidMats,'    KF.s00 = zeros(%.0f,1);\n\n',s.nStateVar);
end

if isfield(s,'KFinit') && isfield(s.KFinit,'Variance')
    fprintf(fidMats,'    sig00 = [...\n');
    for jeq=1:s.nStateVar
        fprintf(fidMats,'       ');
        for jc=1:s.nStateVar
            fprintf(fidMats,' %0.16f',s.KFinit.Variance(jeq,jc));
            if jc==s.nStateVar
                fprintf(fidMats,';\n');
            else
                fprintf(fidMats,',');
            end
        end
    end
    fprintf(fidMats,'    ];\n\n');
    fprintf(fidMats,'    sig00rc = 0;\n');
else
    fprintf(fidMats,'    [sig00,sig00rc] = lyapcsd(REE.G1,REE.G2*REE.G2'');\n');
    fprintf(fidMats,'    sig00 = real(sig00); sig00 = (sig00+sig00'')/2;\n');
    fprintf(fidMats,'    if sig00rc~=0\n');
    fprintf(fidMats,'        if verbose\n');
    fprintf(fidMats,['            fprintf(fid,''Warning: Could not find ',...
                     'unconditional variance!\\n'');\n']);
    fprintf(fidMats,'        end\n\n');
    fprintf(fidMats,'    end\n\n');
end
fprintf(fidMats,'    KF.sig00 = sig00;\n');
fprintf(fidMats,'    KF.sig00rc = sig00rc;\n');
fprintf(fidMats,'    Mats.KF = KF;\n');
fprintf(fidMats,'end\n');


fprintf(fidMats,'\n%% Auxiliary equations matrices\n');
MatNames = {'one','StateVar_t','StateVar_tF','StateVar_tL','ShockVar_t'};
nCols = [1,s.nStateVar,s.nStateVar,s.nStateVar,s.nShockVar];
fprintf(fidMats,'if op.StoreAuxEq\n');
for jM=1:length(MatNames)
    Mj = MatNames{jM};
    SymMats.AuxEq.(Mj) = jacobian(AuxEq,eval(Mj));
    fprintf(fidMats,'    AuxEq.%s = [...\n',Mj);
    for jeq=1:s.nAuxVar
        fprintf(fidMats,'       ');
        for jc=1:nCols(jM)
            fprintf(fidMats,' %s',...
                    char(eval(sprintf('SymMats.AuxEq.%s(jeq,jc)',Mj))));
            if jc==nCols(jM)
                fprintf(fidMats,';\n');
            else
                fprintf(fidMats,',');
            end
        end
    end
    fprintf(fidMats,'        ];\n\n');
end
fprintf(fidMats,'    if op.SolveREE\n');
fprintf(fidMats,['        AuxEq.REE.GBar = ',...
                 'AuxEq.one+AuxEq.StateVar_tF*REE.GBar;\n']);
fprintf(fidMats,['        AuxEq.REE.G0 = ',...
                 'AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1;\n']);
fprintf(fidMats,['        AuxEq.REE.G1 = AuxEq.StateVar_tL;\n']);
fprintf(fidMats,['        AuxEq.REE.G2 = AuxEq.ShockVar_t;\n']);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if op.StoreAuxEq\n');
fprintf(fidMats,'    Mats.AuxEq = AuxEq;\n');
fprintf(fidMats,'end\n');


% close file
fclose(fidMats);

    


%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
