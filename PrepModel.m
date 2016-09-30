function dsge = PrepModel(dsge)

% PrepModel
%
% Analyzes the model and generates code to evaluate the model for a given
% parameter vector.
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge", however in the codes we use that notation, especially in 
%       PrepModel in order to minimize variable and parameter name 
%       incompatibilities.
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

% % check if model already prepared
% if isfield(dsge.Status,Action) && dsge.Statudsge.(Action), return, end

fprintf('\n*** Analyzing DSGE model\n')

% Set Timer
dsge.TimeElapsed.(Action) = toc();

% Check settings
if ~isfield(dsge.Options,'GensysAuthor'), 
    dsge.Options.GensysAuthor = 'CS';
end


%% -------------------------------------------------------------------

%% Prepare variables and equations

fprintf('Generating symbolic variables and systems of equations...\n')

% Parameters
[dsge.n.Param,nc] = size(dsge.Param);
Param.Names = {dsge.Param{:,1}}';
if nc==4
    Param.PrettyNames = Param.Names;
else
    Param.PrettyNames = {dsge.Param{:,5}}';
end
Param.PriorDist = {dsge.Param{:,2}}';
Param.PriorMean = [dsge.Param{:,3}]';
Param.PriorSD = [dsge.Param{:,4}]';
dsge.Param = Param;
vcSym(Param.Names{:})

% NumSolveParam
NumSolveParam = dsge.NumSolveParam;
if isfield(dsge.NumSolveParam,'Names')
    dsge.n.NumSolveParam = length(NumSolveParam.Names);
else
    dsge.n.NumSolveParam = 0;
end
if dsge.n.NumSolveParam>0
    if ~isfield(NumSolveParam,'PrettyNames')
        NumSolveParam.PrettyNames = NumSolveParam.Names;
    end
    if ~isfield(NumSolveParam,'Eq')
        error('NumSolveParam.Eq not defined!')
    end
    nNumSolveEq = length(NumSolveParam.Eq);
    if nNumSolveEq~=dsge.n.NumSolveParam
        error(['Number of param to solve numerically (%.0f) does not match ' ...
               'number of equations (%.0f)!'], dsge.n.NumSolveParam, ...
              nNumSolveEq)
    end
    if ~isfield(NumSolveParam,'Guess')
        fprintf('NumSolveParam.Guess not defined. Ones assumed.')
        NumSolveParam.Guess = ones(dsge.n.NumSolveParam,1);
    end
    vcSym(NumSolveParam.Names{:})
    dsge.NumSolveParam = NumSolveParam;
end

% Auxiliary Parameters
[dsge.n.AuxParam,nc] = size(dsge.AuxParam);
AuxParam.Names = {dsge.AuxParam{:,1}}';
if nc==2
    AuxParam.PrettyNames = AuxParam.Names;
else
    AuxParam.PrettyNames = {dsge.AuxParam{:,3}}';
end
AuxParam.Expressions = {dsge.AuxParam{:,2}}';
dsge.AuxParam = AuxParam;
vcSym(AuxParam.Names{:})

% Observation variables
[dsge.n.ObsVar,nc] = size(dsge.ObsVar);
ObsVar.Names = {dsge.ObsVar{:,1}}';
if nc==1
    ObsVar.PrettyNames = ObsVar.Names;
else
    ObsVar.PrettyNames = {dsge.ObsVar{:,2}}';
end
dsge.ObsVar = ObsVar;
ObsVar_t = sym(zeros(1,dsge.n.ObsVar));
for j=1:dsge.n.ObsVar
    jv = [ObsVar.Names{j},'_t'];
    vcSym(jv)
    ObsVar_t(j) = eval(jv);
end

% State space variables
[dsge.n.StateVar,nc] = size(dsge.StateVar);
StateVar.Names = {dsge.StateVar{:,1}}';
if nc==1
    StateVar.PrettyNames = StateVar.Names;
else
    StateVar.PrettyNames = {dsge.StateVar{:,2}}';
end
dsge.StateVar = StateVar;
StateVar_t = sym(zeros(1,dsge.n.StateVar)); 
StateVar_tF = sym(zeros(1,dsge.n.StateVar)); 
StateVar_tL = sym(zeros(1,dsge.n.StateVar));
for j=1:dsge.n.StateVar
    jv = [StateVar.Names{j},'_t'];
    vcSym(jv,[jv,'F'],[jv,'L'])
    StateVar_t(j) = eval(jv);
    StateVar_tF(j) = eval([jv,'F']);
    StateVar_tL(j) = eval([jv,'L']);
end

% Shocks
[dsge.n.ShockVar,nc] = size(dsge.ShockVar);
ShockVar.Names = {dsge.ShockVar{:,1}}';
if nc==1
    ShockVar.PrettyNames = ShockVar.Names;
else
    ShockVar.PrettyNames = {dsge.ShockVar{:,2}}';
end
dsge.ShockVar = ShockVar;
ShockVar_t = sym(zeros(1,dsge.n.ShockVar));
for j=1:dsge.n.ShockVar
    jv = [ShockVar.Names{j},'_t'];
    vcSym(jv)
    ShockVar_t(j) = eval(jv);
end

% constant variable
vcSym('one')

% Auxiliary variables
[dsge.n.AuxVar,nc] = size(dsge.AuxVar);
AuxVar.Names = {dsge.AuxVar{:,1}}';
if nc==2
    AuxVar.PrettyNames = AuxVar.Names;
else
    AuxVar.PrettyNames = {dsge.AuxVar{:,3}}';
end
dsge.AuxEq = {dsge.AuxVar{:,2}}.';
dsge.AuxVar = AuxVar;
AuxEq = sym(zeros(dsge.n.AuxVar,1));
for j=1:dsge.n.AuxVar
    jv = [AuxVar.Names{j},'_t'];
    eval([jv,' = ',dsge.AuxEq{j},';'])
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
if dsge.n.ObsVar>0
    dsge.n.ObsEq = length(dsge.ObsEq);
    if dsge.n.ObsEq~=dsge.n.ObsVar
        error(['Number of observables (%\.0f) is different from number of ' ...
               'observation equations (%.0f).'],dsge.n.ObsVar,dsge.n.ObsEq)
    end
    ObsEq = sym(zeros(dsge.n.ObsEq,1));
    for j=1:dsge.n.ObsEq
        ObsEq(j) = eval(dsge.ObsEq{j});
    end
end

% Build State equations
dsge.n.StateEq = length(dsge.StateEq);
if dsge.n.StateEq~=dsge.n.StateVar
    error(['Number of state variables (%.0f) is different from number of ' ...
           'state equations (%.0f).'],dsge.n.StateVar,dsge.n.StateEq)
end
StateEq = sym(zeros(dsge.n.StateEq,1));
for j=1:dsge.n.StateEq
    StateEq(j) = eval(dsge.StateEq{j});
end

%% Generate matrices

fprintf('Generating Mats for model evaluation\n')
dsge.FileName.Mats = sprintf('%s_Mats',dsge.Spec);

% Initiate file
fidMats = fopen([dsge.FileName.Mats,'.m'],'wt');
fprintf(fidMats,'function Mats = %s(x,varargin)\n\n',dsge.FileName.Mats);
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

fprintf(fidMats,'\n%% Initiate Status\n');
fprintf(fidMats,'Mats.Status = 1;\n');
fprintf(fidMats,'Mats.StatusMessage = '''';\n');

if dsge.n.Param>0
    fprintf(fidMats,'\n%% Map parameters\n');
    for j=1:dsge.n.Param
        fprintf(fidMats,'%s = x(%.0f);\n',Param.Names{j},j);
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:dsge.n.Param
        fprintf(fidMats,'    Mats.Param.%s = x(%.0f);\n',Param.Names{j},j);
    end
    fprintf(fidMats,'end\n');
end

if dsge.n.AuxParam>0 || dsge.n.NumSolveParam>0
    fprintf(fidMats,'\n%% Initialize auxiliary parameters\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'%s = [];\n',NumSolveParam.Names{j});
    end
    for j=1:dsge.n.AuxParam
        fprintf(fidMats,'%s = [];\n',AuxParam.Names{j});
    end
end

if dsge.n.NumSolveParam>0
    fprintf(fidMats,'\n%% NumSolve parameters\n');
    fprintf(fidMats,'function f=NumSolveEq(x)\n');
    fprintf(fidMats,'    for jx=1:size(x,2)\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'        %s = x(%.0f,jx);\n',NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'        EvalAuxParam\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'        f(%.0f,jx) = %s;\n',j,NumSolveParam.Eq{j});
    end
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'NumSolveGuess(%.0f,1) = %.16f;\n',...
                j,NumSolveParam.Guess(j));
    end
    fprintf(fidMats,['[NumSolveSolution,NumSolveRC] = csolvevb(@NumSolveEq,' ...
                     'NumSolveGuess,[],1e-10,1000);']);
    fprintf(fidMats,'Mats.NumSolveParamRC = NumSolveRC;\n');
    fprintf(fidMats,'if NumSolveRC~=0 && verbose\n');
    fprintf(fidMats,'    Mats.Status = 0;\n');
    txt = 'NumSolveParam solution not normal!';
    fprintf(fidMats,'    Mats.StatusMessage = [Mats.StatusMessage,''%s\\n''];\n',...
            txt);
    fprintf(fidMats,'    fprintf(fid,''Warning: %s\\n'');\n',txt);
    fprintf(fidMats,'end\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'%s = NumSolveSolution(%.0f);\n',...
                NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:dsge.n.NumSolveParam
        fprintf(fidMats,'    Mats.NumSolveParam.%1$s = %1$s;\n',...
                NumSolveParam.Names{j});
    end
    fprintf(fidMats,'end\n');
end

if dsge.n.AuxParam>0
    fprintf(fidMats,'\n%% Map auxiliary parameters\n');
    if dsge.n.NumSolveParam>0
        fprintf(fidMats,'function EvalAuxParam \n');
    end
    for j=1:dsge.n.AuxParam
        fprintf(fidMats,'%s = %s;\n',AuxParam.Names{j},AuxParam.Expressions{j});
    end
    if dsge.n.NumSolveParam>0
        fprintf(fidMats,'end \n');
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:dsge.n.AuxParam
        fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',AuxParam.Names{j});
    end
    fprintf(fidMats,'end\n');
end

if dsge.n.ObsVar>0
    fprintf(fidMats,'\n%% Observation equations\n');
    H0 = -jacobian(ObsEq,ObsVar_t);
    SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
    SymMats.ObsEq.H = H0\jacobian(ObsEq,StateVar_t);
    MatNames = fieldnames(SymMats.ObsEq);
    nCols = [1,dsge.n.StateVar];
    fprintf(fidMats,'if op.StoreObsEq || op.StoreKF\n');
    for jM=1:length(MatNames)
        fprintf(fidMats,'    ObsEq.%s = [...\n',MatNames{jM});
        for jeq=1:dsge.n.ObsVar
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
end

fprintf(fidMats,'\n%% State equation matrices\n');
SymMats.StateEq.GammaBar = jacobian(StateEq,one);
SymMats.StateEq.Gamma0 = -jacobian(StateEq,StateVar_tF);
SymMats.StateEq.Gamma1 = jacobian(StateEq,StateVar_t);
SymMats.StateEq.Gamma4 = jacobian(StateEq,StateVar_tL);
SymMats.StateEq.Gamma2 = jacobian(StateEq,ShockVar_t);
MatNames = fieldnames(SymMats.StateEq);
nCols = [1,dsge.n.StateVar,dsge.n.StateVar,dsge.n.StateVar,dsge.n.ShockVar];
for jM=1:length(MatNames)
    fprintf(fidMats,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:dsge.n.StateVar
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
fprintf(fidMats,'StateEq.Gamma3 = eye(%.0f);\n\n',dsge.n.StateVar);
fprintf(fidMats,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',...
        dsge.n.StateVar);
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
        dsge.Options.GensysAuthor);
fprintf(fidMats,'    Mats.REE = REE;\n');
fprintf(fidMats,'    if ~all(REE.eu==1);\n');
fprintf(fidMats,'        Mats.Status = 0;\n');
fprintf(fidMats,['        Mats.StatusMessage = [Mats.StatusMessage,''REE ' ...
                 'solution not normal!\\n''];\n']);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'end\n');

if dsge.n.ObsVar>0
    fprintf(fidMats,'\n%% Kalman Filter matrices\n');
    fprintf(fidMats,'if op.StoreKF\n');
    fprintf(fidMats,'    if all(Mats.REE.GBar(:)==0)\n');
    fprintf(fidMats,'        KF.StateVarBar = zeros(%.0f,1);\n',...
            dsge.n.StateVar);
    fprintf(fidMats,'    else\n');
    fprintf(fidMats,...
            '        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
            dsge.n.StateVar);
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,...
            '    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

    if isfield(dsge,'KFinit') && isfield(dsge.KFinit,'State')
        fprintf(fidMats,'    s00 = [...\n');
        for jeq=1:dsge.n.StateVar
            fprintf(fidMats,'        %.16f;\n',dsge.KFinit.State(jeq));
        end
        fprintf(fidMats,'        ];\n\n');
    else
        fprintf(fidMats,'    KF.s00 = zeros(%.0f,1);\n\n',dsge.n.StateVar);
    end

    if isfield(dsge,'KFinit') && isfield(dsge.KFinit,'Variance')
        fprintf(fidMats,'    sig00 = [...\n');
        for jeq=1:dsge.n.StateVar
            fprintf(fidMats,'       ');
            for jc=1:dsge.n.StateVar
                fprintf(fidMats,' %0.16f',dsge.KFinit.Variance(jeq,jc));
                if jc==dsge.n.StateVar
                    fprintf(fidMats,';\n');
                else
                    fprintf(fidMats,',');
                end
            end
        end
        fprintf(fidMats,'    ];\n\n');
        fprintf(fidMats,'    sig00rc = 0;\n');
    else
        fprintf(fidMats,...
                '    [sig00,sig00rc] = lyapcsd(REE.G1,REE.G2*REE.G2'');\n');
        fprintf(fidMats,...
                '    sig00 = real(sig00); sig00 = (sig00+sig00'')/2;\n');
        fprintf(fidMats,'    if sig00rc~=0\n');
        txt = 'Could not find unconditional variance!';
        fprintf(fidMats,...
                '        Mats.StatusMessage = [Mats.StatusMessage,''%s\\n''];\n',...
                txt);
        fprintf(fidMats,'        if verbose\n');
        fprintf(fidMats,...
                '            fprintf(fid,''Warning: %s\\n'');\n',txt);
        fprintf(fidMats,'        end\n\n');
        fprintf(fidMats,'    end\n\n');
    end
    fprintf(fidMats,'    KF.sig00 = sig00;\n');
    fprintf(fidMats,'    KF.sig00rc = sig00rc;\n');
    fprintf(fidMats,'    Mats.KF = KF;\n');
    fprintf(fidMats,'end\n');
end

if dsge.n.AuxVar>0
    fprintf(fidMats,'\n%% Auxiliary equations matrices\n');
    MatNames = {'one','StateVar_t','StateVar_tF','StateVar_tL','ShockVar_t'};
    nCols = [1,dsge.n.StateVar,dsge.n.StateVar,dsge.n.StateVar,dsge.n.ShockVar];
    fprintf(fidMats,'if op.StoreAuxEq || op.StoreAuxREE\n');
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        SymMats.AuxEq.(Mj) = jacobian(AuxEq,eval(Mj));
        fprintf(fidMats,'    AuxEq.%s = [...\n',Mj);
        for jeq=1:dsge.n.AuxVar
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
    fprintf(fidMats,'        if ~isempty(REE.G1)\n');
    fprintf(fidMats,['            AuxREE.GBar = ',...
                     'AuxEq.one+AuxEq.StateVar_tF*REE.GBar',...
                     '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.GBar',...
                     ';\n']);
    fprintf(fidMats,['            AuxREE.G1 = AuxEq.StateVar_tL',...
                     '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.G1;\n']);
    fprintf(fidMats,['            AuxREE.G2 = AuxEq.ShockVar_t',...
                     '+(AuxEq.StateVar_t+AuxEq.StateVar_tF*REE.G1)*REE.G2;\n']);
    fprintf(fidMats,'        else\n');
    fprintf(fidMats,'            AuxREE.GBar = [];\n');
    fprintf(fidMats,'            AuxREE.G1 = [];\n');
    fprintf(fidMats,'            AuxREE.G2 = [];\n');
    fprintf(fidMats,'        end\n');
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    fprintf(fidMats,'if op.StoreAuxEq\n');
    fprintf(fidMats,'    Mats.AuxEq = AuxEq;\n');
    fprintf(fidMats,'end\n');
    fprintf(fidMats,'if op.StoreAuxREE\n');
    fprintf(fidMats,'    Mats.AuxREE = AuxREE;\n');
    fprintf(fidMats,'end\n');
end

% close file
fprintf(fidMats,'end\n');
fclose(fidMats);

    


%% -------------------------------------------------------------------

%% Finish up
dsge.Status.(Action) = 1;
dsge.TimeElapsed.(Action) = toc-dsge.TimeElapsed.(Action);
fprintf('\n%s %s\n\n',Action,vctoc([],dsge.TimeElapsed.(Action)))

%% -------------------------------------------------------------------
