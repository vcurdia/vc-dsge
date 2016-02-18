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
if ~isfield(s.Options,'GensysAuthor'), s.Options.GensysAuthor = 'CS';end


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
n.Param = length(Param.Names);
for j=1:n.Param
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
n.AuxParam = length(AuxParam.Names);
for j=1:n.AuxParam
    eval(['syms ',AuxParam.Names{j}]);
end

% Observation variables
n.ObsVar = length(s.ObsVar);
ObsVar_t = sym(zeros(1,n.ObsVar));
for j=1:n.ObsVar
    jv = [s.ObsVar{j},'_t'];
    eval(['syms ',jv]);
    ObsVar_t(j) = eval(jv);
end

% State space variables
n.StateVar = length(s.StateVar);
StateVar_t = sym(zeros(1,n.StateVar)); 
StateVar_tF = sym(zeros(1,n.StateVar)); 
StateVar_tL = sym(zeros(1,n.StateVar));
for j=1:n.StateVar
    jv = [s.StateVar{j},'_t'];
    eval(sprintf('syms %1$s %1$sF %1$sL',jv))
    StateVar_t(j) = eval(jv);
    StateVar_tF(j) = eval([jv,'F']);
    StateVar_tL(j) = eval([jv,'L']);
end

% Shocks
n.ShockVar = length(s.ShockVar);
ShockVar_t = sym(zeros(1,n.ShockVar));
for j=1:n.ShockVar
    jv = [s.ShockVar{j},'_t'];
    eval(['syms ',jv]);
    ShockVar_t(j) = eval(jv);
end

% constant variable
syms one

% Auxiliary variables
s.AuxEq = {s.AuxVar{:,2}}.';
s.AuxVar = {s.AuxVar{:,1}}.';
n.AuxVar = length(s.AuxVar);
AuxEq = sym(zeros(n.AuxVar,1));
for j=1:n.AuxVar
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
n.ObsEq = length(s.ObsEq);
if n.ObsEq~=n.ObsVar
    error(['Number of observables (%\.0f) is different from number of ' ...
           'observation equations (%.0f).'],n.ObsVar,n.ObsEq)
end
ObsEq = sym(zeros(n.ObsEq,1));
for j=1:n.ObsEq
    ObsEq(j) = sym(s.ObsEq{j});
end

% Build State equations
n.StateEq = length(s.StateEq);
if n.StateEq~=n.StateVar
    error(['Number of state variables (%\.0f) is different from number of ' ...
           'state equations (%.0f).'],n.StateVar,n.StateEq)
end
StateEq = sym(zeros(n.StateEq,1));
for j=1:n.StateEq
    StateEq(j) = sym(s.StateEq{j});
end

s.n = n;


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
fprintf(fidMats,'op.StoreAuxREE = 1;\n');
fprintf(fidMats,'op.SolveREE = 1;\n');
fprintf(fidMats,'op.fid = 1;\n');
fprintf(fidMats,'op.verbose = 0;\n');
fprintf(fidMats,'op.gensys = {};\n');

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
for j=1:n.Param
    fprintf(fidMats,'%s = x(%.0f);\n',Param.Names{j},j);
end
fprintf(fidMats,'if op.StoreParam\n');
for j=1:n.Param
    fprintf(fidMats,'    Mats.Param.%s = x(%.0f);\n',Param.Names{j},j);
end
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Map auxiliary parameters\n');
for j=1:n.AuxParam
    fprintf(fidMats,'%s = %s;\n',AuxParam.Names{j},AuxParam.Expressions{j});
end
fprintf(fidMats,'if op.StoreParam\n');
for j=1:n.AuxParam
    fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',AuxParam.Names{j});
end
fprintf(fidMats,'end\n');

fprintf(fidMats,'\n%% Observation equations\n');
H0 = -jacobian(ObsEq,ObsVar_t);
SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
SymMats.ObsEq.H = H0\jacobian(ObsEq,StateVar_t);
MatNames = fieldnames(SymMats.ObsEq);
nCols = [1,n.StateVar];
fprintf(fidMats,'if op.StoreObsEq || op.StoreKF\n');
for jM=1:length(MatNames)
    fprintf(fidMats,'    ObsEq.%s = [...\n',MatNames{jM});
    for jeq=1:n.ObsVar
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
nCols = [1,n.StateVar,n.StateVar,n.StateVar,n.ShockVar];
for jM=1:length(MatNames)
    fprintf(fidMats,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:n.StateVar
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
fprintf(fidMats,'StateEq.Gamma3 = eye(%.0f);\n\n',n.StateVar);
fprintf(fidMats,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',n.StateVar);
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
        '        ''%s'',op.fid,op.verbose,op.gensys{:});\n',...
        s.Options.GensysAuthor);
fprintf(fidMats,'    Mats.REE = REE;\n');
fprintf(fidMats,'end\n');


fprintf(fidMats,'\n%% Kalman Filter matrices\n');
fprintf(fidMats,'if op.StoreKF\n');
fprintf(fidMats,'    if all(Mats.REE.GBar(:)==0)\n');
fprintf(fidMats,'        KF.StateVarBar = zeros(%.0f,1);\n',n.StateVar);
fprintf(fidMats,'    else\n');
fprintf(fidMats,'        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
        n.StateVar);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

if isfield(s,'KFinit') && isfield(s.KFinit,'State')
    fprintf(fidMats,'    s00 = [...\n');
    for jeq=1:n.StateVar
        fprintf(fidMats,'        %.16f;\n',s.KFinit.State(jeq));
    end
    fprintf(fidMats,'        ];\n\n');
else
    fprintf(fidMats,'    KF.s00 = zeros(%.0f,1);\n\n',n.StateVar);
end

if isfield(s,'KFinit') && isfield(s.KFinit,'Variance')
    fprintf(fidMats,'    sig00 = [...\n');
    for jeq=1:n.StateVar
        fprintf(fidMats,'       ');
        for jc=1:n.StateVar
            fprintf(fidMats,' %0.16f',s.KFinit.Variance(jeq,jc));
            if jc==n.StateVar
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
nCols = [1,n.StateVar,n.StateVar,n.StateVar,n.ShockVar];
fprintf(fidMats,'if op.StoreAuxEq || op.StoreAuxREE\n');
for jM=1:length(MatNames)
    Mj = MatNames{jM};
    SymMats.AuxEq.(Mj) = jacobian(AuxEq,eval(Mj));
    fprintf(fidMats,'    AuxEq.%s = [...\n',Mj);
    for jeq=1:n.AuxVar
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
fprintf(fidMats,'    if op.SolveREE && op.StoreAuxREE\n');
fprintf(fidMats,['        AuxREE.GBar = ',...
                 'AuxEq.one+AuxEq.StateVar_tF*REE.GBar',...
                 '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.GBar;\n']);
fprintf(fidMats,['        AuxREE.G1 = AuxEq.StateVar_tL',...
                 '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.G1;\n']);
fprintf(fidMats,['        AuxREE.G2 = AuxEq.ShockVar_t',...
                 '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.G2;\n']);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if op.StoreAuxEq\n');
fprintf(fidMats,'    Mats.AuxEq = AuxEq;\n');
fprintf(fidMats,'end\n');
fprintf(fidMats,'if op.StoreAuxREE\n');
fprintf(fidMats,'    Mats.AuxREE = AuxREE;\n');
fprintf(fidMats,'end\n');


% close file
fclose(fidMats);

    


%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------
