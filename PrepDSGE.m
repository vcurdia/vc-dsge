function obj = PrepDSGE(obj)

% PrepDSGE
%
% Analyzes the DSGE structure and generates code to evaluate the DSGE for a 
% given parameter vector.
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "obj", however in the codes we use that notation, especially in 
%       PrepDSGE in order to minimize variable and parameter name 
%       incompatibilities.
%
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

Action = 'PrepDSGE';

% % check if model already prepared
% if isfield(obj.Status,Action) && obj.Statudsge.(Action), return, end

fprintf('\n*** Analyzing DSGE model\n')

% Set Timer
obj.TimeElapsed.(Action) = toc();

% Check settings
if ~isfield(obj.Options,'GensysAuthor'), 
    obj.Options.GensysAuthor = 'CS';
end


%% -------------------------------------------------------------------

%% Prepare parameters, variables and equations

fprintf('Generating symbolic variables and systems of equations...\n')

%% Parameters
pp = obj.Param;
pList = {'Est','Set','NumSolve','Aux'};
p = struct;
for j=1:length(pList)
    jp = pList{j};
    if isfield(obj.Param,jp)
        [obj.n.Param.(jp),nc] = size(pp.(jp));
    else
        obj.n.Param.(jp) = 0;
    end
    obj.Exist.Param.(jp) = obj.n.Param.(jp)>0;
    if obj.Exist.Param.(jp)
        p.(jp) = struct;
        p.(jp).Names = {pp.(jp){:,1}}';
        if nc==3+2*strcmp(jp,'Est')
            p.(jp).PrettyNames = {pp.(jp){:,nc}}';
        else
            p.(jp).PrettyNames = p.(jp).Names;
        end
        vcSym(p.(jp).Names{:})
    end
end

if obj.Exist.Param.Est
    p.Est.PriorDist = {pp.Est{:,2}}';
    p.Est.PriorMean = [pp.Est{:,3}]';
    p.Est.PriorSD = [pp.Est{:,4}]';
end

if obj.Exist.Param.Set
    p.Set.Values = [pp.Set{:,2}]';
end

if obj.Exist.Param.NumSolve
    p.NumSolve.Guess = [pp.NumSolve{:,2}];
    if ~isfield(pp,'NumSolveEq')
        error('NumSolveEq not defined!')
    else
        p.NumSolve.Eq = pp.NumSolveEq;
    end
    nNumSolveEq = length(pp.NumSolveEq);
    if nNumSolveEq~=obj.n.Param.NumSolve
        error(['Number of param to solve numerically (%.0f) does not match ' ...
               'number of equations (%.0f)!'], obj.n.Param.NumSolve, ...
              nNumSolveEq)
    end
end

if obj.Exist.Param.Aux
    p.Aux.Expressions = {pp.Aux{:,2}}';
end

obj.Param = p;


%% Variables

vv = obj.Var;
vList = {'Obs','State','Shock','Aux'};
v = struct;
vcSym('one')
Var = struct;
for j=1:length(vList)
    jv = vList{j};
    if isfield(obj.Var,jv)
        [obj.n.Var.(jv),nc] = size(vv.(jv));
    else
        obj.n.Var.(jv) = 0;
    end
    obj.Exist.Var.(jv) = obj.n.Var.(jv)>0;
    if obj.Exist.Var.(jv)
        v.(jv) = struct;
        v.(jv).Names = {vv.(jv){:,1}}';
        if nc==2+strcmp(jp,'Aux')
            v.(jv).PrettyNames = {vv.(jv){:,nc}}';
        else
            v.(jv).PrettyNames = v.(jv).Names;
        end
        Var.(jv).t = sym(zeros(1,obj.n.Var.(jv)));
        for jj=1:obj.n.Var.(jv)
            vj = [v.(jv).Names{jj},'_t'];
            vcSym(vj)
            Var.(jv).t(j) = eval(vj);
        end
    end
end

if obj.Exist.Var.State
    Var.State.tF = sym(zeros(1,obj.n.Var.State)); 
    Var.State.tL = sym(zeros(1,obj.n.Var.State));
    for j=1:obj.n.Var.State
        jv = [v.State.Names{j},'_t'];
        vcSym([jv,'F'],[jv,'L'])
        Var.State.tF(j) = eval([jv,'F']);
        Var.State.tL(j) = eval([jv,'L']);
    end
end

if obj.Exist.Var.Aux
    obj.n.Eq.Aux = obj.n.Var.Aux;
    obj.Eq.Aux = {vv.Aux{:,2}}.';
    EqAux = sym(zeros(obj.n.Var.Aux,1));
    for j=1:obj.n.Eq.Aux
        jv = [v.Aux.Names{j},'_t'];
        eval([jv,' = ',obj.Eq.Aux{j},';'])
        EqAux(j) = eval(jv);
% If the expression has no leads then can define a lead for it
        if all(jacobian(EqAux(j),Var.State.tF)==0)
            eval([jv,'F = subs(',jv,',[Var.State.t,Var.State.tL],',...
                  '[Var.State.tF,Var.State.t]);'])
% If expreassion has no leads or lags then can define lag for it. 
% Notice that it does not make sense to define a lag if there are leads in it,
% and that's why the check for lags is inside the check for leads
            if all(jacobian(EqAux(j),Var.State.tL)==0)
                eval([jv,'L = subs(',jv,',[Var.State.tF,Var.State.t],',...
                      '[Var.State.t,Var.State.tL]);'])
            end
        end
    end
end

% Build Observation equations
if obj.n.Var.Obs>0
    obj.n.Eq.Obs = length(obj.Eq.Obs);
    if obj.n.Eq.Obs~=obj.n.Var.Obs
        error(['Number of observables (%\.0f) is different from number of ' ...
               'observation equations (%.0f).'],obj.n.Var.Obs,obj.n.Eq.Obs)
    end
    EqObs = sym(zeros(obj.n.Eq.Obs,1));
    for j=1:obj.n.Eq.Obs
        EqObs(j) = eval(obj.Eq.Obs{j});
    end
end

% Build State equations
obj.n.Eq.State = length(obj.Eq.State);
if obj.n.Eq.State~=obj.n.Var.State
    error(['Number of state variables (%.0f) is different from number of ' ...
           'state equations (%.0f).'],obj.n.Var.State,obj.n.Eq.State)
end
EqState = sym(zeros(obj.n.Eq.State,1));
for j=1:obj.n.Eq.State
    EqState(j) = eval(obj.Eq.State{j});
end

keyboard


%% Generate matrices

fprintf('Generating Mats for model evaluation\n')
obj.FileName.Mats = sprintf('%s_Mats',obj.Spec);

% Initiate file
fidMats = fopen([obj.FileName.Mats,'.m'],'wt');
fprintf(fidMats,'function Mats = %s(x,varargin)\n\n',obj.FileName.Mats);
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

if obj.n.Param>0
    fprintf(fidMats,'\n%% Map parameters\n');
    for j=1:obj.n.Param
        fprintf(fidMats,'%s = x(%.0f);\n',Param.Names{j},j);
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:obj.n.Param
        fprintf(fidMats,'    Mats.Param.%s = x(%.0f);\n',Param.Names{j},j);
    end
    fprintf(fidMats,'end\n');
end

if obj.n.AuxParam>0 || obj.n.NumSolveParam>0
    fprintf(fidMats,'\n%% Initialize auxiliary parameters\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'%s = [];\n',NumSolveParam.Names{j});
    end
    for j=1:obj.n.AuxParam
        fprintf(fidMats,'%s = [];\n',AuxParam.Names{j});
    end
end

if obj.n.NumSolveParam>0
    fprintf(fidMats,'\n%% NumSolve parameters\n');
    fprintf(fidMats,'function f=NumSolveEq(x)\n');
    fprintf(fidMats,'    for jx=1:size(x,2)\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'        %s = x(%.0f,jx);\n',NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'        EvalAuxParam\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'        f(%.0f,jx) = %s;\n',j,NumSolveParam.Eq{j});
    end
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'NumSolveGuess(%.0f,1) = %.16f;\n',...
                j,NumSolveParam.Guess(j));
    end
    fprintf(fidMats,['[NumSolveSolution,NumSolveRC] = csolvevb(@NumSolveEq,' ...
                     'NumSolveGuess,[],1e-10,1000);']);
    fprintf(fidMats,'Mats.NumSolveParamRC = NumSolveRC;\n');
    fprintf(fidMats,'if NumSolveRC~=0\n');
    fprintf(fidMats,'    Mats.Status = 0;\n');
    txt = 'NumSolveParam solution not normal!';
    fprintf(fidMats,'    Mats.StatusMessage = [Mats.StatusMessage,''%s\\n''];\n',...
            txt);
    fprintf(fidMats,'    if verbose\n');
    fprintf(fidMats,'        fprintf(fid,''Warning: %s\\n'');\n',txt);
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'%s = NumSolveSolution(%.0f);\n',...
                NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:obj.n.NumSolveParam
        fprintf(fidMats,'    Mats.NumSolveParam.%1$s = %1$s;\n',...
                NumSolveParam.Names{j});
    end
    fprintf(fidMats,'end\n');
end

if obj.n.AuxParam>0
    fprintf(fidMats,'\n%% Map auxiliary parameters\n');
    if obj.n.NumSolveParam>0
        fprintf(fidMats,'function EvalAuxParam \n');
    end
    for j=1:obj.n.AuxParam
        fprintf(fidMats,'%s = %s;\n',AuxParam.Names{j},AuxParam.Expressions{j});
    end
    if obj.n.NumSolveParam>0
        fprintf(fidMats,'end \n');
    end
end

% Incorporate NumSolveParam into AuxParam
if obj.n.NumSolveParam>0
    obj.n.AuxParam = obj.n.AuxParam + obj.n.NumSolveParam;
    obj.AuxParam.Names = [obj.AuxParam.Names;obj.NumSolveParam.Names];
    obj.AuxParam.PrettyNames = ...
        [obj.AuxParam.PrettyNames;obj.NumSolveParam.PrettyNames];
    obj.AuxParam.Expressions(obj.n.AuxParam+1-(1:obj.n.NumSolveParam)) = ...
                        {'NaN'};
end
if obj.n.AuxParam>0
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:obj.n.AuxParam
        fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',...
                obj.AuxParam.Names{j});
    end
    fprintf(fidMats,'end\n');
end

if obj.n.ObsVar>0
    fprintf(fidMats,'\n%% Observation equations\n');
    H0 = -jacobian(ObsEq,ObsVar_t);
    SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
    SymMats.ObsEq.H = H0\jacobian(ObsEq,StateVar_t);
    MatNames = fieldnames(SymMats.ObsEq);
    nCols = [1,obj.n.StateVar];
    fprintf(fidMats,'if op.StoreObsEq || op.StoreKF\n');
    for jM=1:length(MatNames)
        fprintf(fidMats,'    ObsEq.%s = [...\n',MatNames{jM});
        for jeq=1:obj.n.ObsVar
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
nCols = [1,obj.n.StateVar,obj.n.StateVar,obj.n.StateVar,obj.n.ShockVar];
for jM=1:length(MatNames)
    fprintf(fidMats,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:obj.n.StateVar
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
fprintf(fidMats,'StateEq.Gamma3 = eye(%.0f);\n\n',obj.n.StateVar);
fprintf(fidMats,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',...
        obj.n.StateVar);
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
        obj.Options.GensysAuthor);
fprintf(fidMats,'    Mats.REE = REE;\n');
fprintf(fidMats,'    if ~all(REE.eu==1);\n');
fprintf(fidMats,'        Mats.Status = 0;\n');
fprintf(fidMats,['        Mats.StatusMessage = [Mats.StatusMessage,''REE ' ...
                 'solution not normal!\\n''];\n']);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'end\n');

if obj.n.ObsVar>0
    fprintf(fidMats,'\n%% Kalman Filter matrices\n');
    fprintf(fidMats,'if op.StoreKF\n');
    fprintf(fidMats,'    if all(Mats.REE.GBar(:)==0)\n');
    fprintf(fidMats,'        KF.StateVarBar = zeros(%.0f,1);\n',...
            obj.n.StateVar);
    fprintf(fidMats,'    else\n');
    fprintf(fidMats,...
            '        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
            obj.n.StateVar);
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,...
            '    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

    if isfield(obj,'KFinit') && isfield(obj.KFinit,'State')
        fprintf(fidMats,'    s00 = [...\n');
        for jeq=1:obj.n.StateVar
            fprintf(fidMats,'        %.16f;\n',obj.KFinit.State(jeq));
        end
        fprintf(fidMats,'        ];\n\n');
    else
        fprintf(fidMats,'    KF.s00 = zeros(%.0f,1);\n\n',obj.n.StateVar);
    end

    if isfield(obj,'KFinit') && isfield(obj.KFinit,'Variance')
        fprintf(fidMats,'    sig00 = [...\n');
        for jeq=1:obj.n.StateVar
            fprintf(fidMats,'       ');
            for jc=1:obj.n.StateVar
                fprintf(fidMats,' %0.16f',obj.KFinit.Variance(jeq,jc));
                if jc==obj.n.StateVar
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

if obj.n.AuxVar>0
    fprintf(fidMats,'\n%% Auxiliary equations matrices\n');
    MatNames = {'one','StateVar_t','StateVar_tF','StateVar_tL','ShockVar_t'};
    nCols = [1,obj.n.StateVar,obj.n.StateVar,obj.n.StateVar,obj.n.ShockVar];
    fprintf(fidMats,'if op.StoreAuxEq || op.StoreAuxREE\n');
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        SymMats.AuxEq.(Mj) = jacobian(AuxEq,eval(Mj));
        fprintf(fidMats,'    AuxEq.%s = [...\n',Mj);
        for jeq=1:obj.n.AuxVar
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
obj.TimeElapsed.(Action) = toc-obj.TimeElapsed.(Action);
fprintf('\n%s %s\n\n',Action,vctoc([],obj.TimeElapsed.(Action)))

%% -------------------------------------------------------------------
