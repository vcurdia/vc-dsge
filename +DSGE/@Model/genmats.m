function genmats(obj)

% genmats
%
% Analyzes the DSGE structure and generates code to evaluate the DSGE for a 
% given parameter vector.
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
%
% See also:
% DSGE, SetupMyDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia

%% Preamble
fprintf('\n*** Generating DSGE mats\n')
ttName = 'GenMats';
obj.TimeTracker.start(ttName)

%% basic check
if (obj.StateVar.N==0) || isempty(obj.StateEq)
    error('Cannot proceed without specifying state variables and equations')
end

%% Sym Params
list = {'','NumSolve','Compound'};
% sList = struct;
for j=1:length(list)
    jstr = [list{j},'Param'];
    if obj.(jstr).N>0, vcsym(obj.(jstr).Names{:}), end
%     sList.(jstr) = sym(zeros(1,obj.(jstr).N));
%     for jp=1:obj.(jstr).N
%         sList.(jstr)(jp) = eval(obj.(jstr).Names{jp});
%     end
end
% SymParam = sList.Param;
% SymNumsolveParam = sList.NumSolveParam;
% SymCompoundParam = sList.CompoundParam;
% SymParamAll = [SymParam,SymNumsolveParam,SymCompoundParam];

% %% Constant
% vcsym('one')

%% Obs Var
ObsVar_t = sym(zeros(1,obj.ObsVar.N)); 
for j=1:obj.ObsVar.N
    vj = [obj.ObsVar.Names{j},'_t'];
    vcsym(vj)
    ObsVar_t(j) = eval(vj);
end

%% State Var
StateVar_t = sym(zeros(1,obj.StateVar.N)); 
StateVar_tF = sym(zeros(1,obj.StateVar.N)); 
StateVar_tL = sym(zeros(1,obj.StateVar.N)); 
for j=1:obj.StateVar.N
    vj = [obj.StateVar.Names{j},'_t'];
    vcsym(vj,[vj,'F'],[vj,'L'])
    StateVar_t(j) = eval(vj);
    StateVar_tF(j) = eval([vj,'F']);
    StateVar_tL(j) = eval([vj,'L']);
end

%% Shock Var
ShockVar_t = sym(zeros(1,obj.ShockVar.N)); 
for j=1:obj.ShockVar.N
    vj = [obj.ShockVar.Names{j},'_t'];
    vcsym(vj)
    ShockVar_t(j) = eval(vj);
end

%% Aux Var and Eq
nV = obj.AuxVar.N;
if nV>0
    AuxEq = sym(zeros(nV,1));
    for j=1:nV
        vj = [obj.AuxVar.Names{j},'_t'];
        eval([vj,' = ',obj.AuxEq{j},';'])
        AuxEq(j) = eval(vj);
% If the expression has no leads then can define a lead for it
        if all(jacobian(AuxEq(j),StateVar_tF)==0)
            eval([vj,'F = subs(',vj,',[StateVar_t,StateVar_tL],',...
                  '[StateVar_tF,StateVar_t]);'])
% If expression has no leads or lags then can define lag for it. 
% Notice that it does not make sense to define a lag if there are leads in it,
% and that's why the check for lags is inside the check for leads
            if all(jacobian(AuxEq(j),StateVar_tL)==0)
                eval([vj,'L = subs(',vj,',[StateVar_tF,StateVar_t],',...
                      '[StateVar_t,StateVar_tL]);'])
            end
        end
    end
end

%% Build Observation equations
if obj.ObsVar.N>0
    ObsEq = sym(zeros(obj.ObsVar.N,1));
    for j=1:obj.ObsVar.N
        ObsEq(j) = eval(obj.ObsEq{j});
    end
end

%% Build State equations
StateEq = sym(zeros(obj.StateVar.N,1));
for j=1:obj.StateVar.N
    StateEq(j) = eval(obj.StateEq{j});
end


%% Generate matrices

fprintf('Generating code to evaluate model Mats\n')

MatsFN = sprintf('%s_Mats',obj.Name);

% Initiate file
fid = fopen([MatsFN,'.m'],'wt');
fprintf(fid,'function Mats = %s(x,varargin)\n\n',MatsFN);
fprintf(fid,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs\n',clock);

fprintf(fid,'\n%% Default options\n');
fprintf(fid,'op.StoreParam = 1;\n');
fprintf(fid,'op.StoreStateEq = 1;\n');
fprintf(fid,'op.StoreObsEq = 1;\n');
fprintf(fid,'op.StoreKF = 1;\n');
fprintf(fid,'op.StoreAuxEq = 1;\n');
fprintf(fid,'op.StoreAuxREE = 1;\n');
fprintf(fid,'op.SolveREE = 1;\n');
fprintf(fid,'op.fid = 1;\n');
fprintf(fid,'op.verbose = 0;\n');
fprintf(fid,'op.GensysAuthor = ''%s'';\n',obj.GensysAuthor);
fprintf(fid,'op.gensys = {};\n');

fprintf(fid,'\n%% Update options\n');
fprintf(fid,'op = updateoptions(op,varargin{:});\n');

fprintf(fid,'\n%% Initiate Status\n');
fprintf(fid,'Mats.Status = 1;\n');
fprintf(fid,'Mats.StatusMessage = '''';\n');

fprintf(fid,'\n%% Map parameters\n');
for j=1:obj.Param.N
    fprintf(fid,'%s = x(%.0f);\n',obj.Param.Names{j},j);
end
fprintf(fid,'if op.StoreParam\n');
fprintf(fid,'    Mats.Param = x;\n');
% for j=1:obj.Param.N
%     fprintf(fid,'    Mats.Param.%1$s = %1$s;\n',obj.Param.Names{j});
% end
fprintf(fid,'end\n');

if obj.CompoundParam.N>0 || obj.NumSolveParam.N>0
    fprintf(fid,'\n%% Initialize compound parameters\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'%s = [];\n',obj.NumSolveParam.Names{j});
    end
    for j=1:obj.CompoundParam.N
        fprintf(fid,'%s = [];\n',obj.CompoundParam.Names{j});
    end
end

if obj.NumSolveParam.N>0
    fprintf(fid,'\n%% NumSolve parameters\n');
    fprintf(fid,'function f=NumSolveEq(x)\n');
    fprintf(fid,'    for jx=1:size(x,2)\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'        %s = x(%.0f,jx);\n',...
                obj.NumSolveParam.Names{j},j);
    end
    fprintf(fid,'        EvalCompoundParam\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'        f(%.0f,jx) = %s;\n',j,...
                obj.NumSolveEq{j});
    end
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    fprintf(fid,'NumSolveGuess = [...\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'    %.16f;\n',obj.NumSolveParam.Values(j));
    end
    fprintf(fid,'    ];\n');
%     fprintf(fid,['[NumSolveSolution,NumSolveRC] = csolvevb(@NumSolveEq,' ...
%                      'NumSolveGuess,[],%e,%.0f);\n'],...
%             obj.NumPrecision,obj.NumSolveMaxIterations);
    fprintf(fid,['NumSolveOptions = optimoptions(@fsolve);\n']);
    fprintf(fid,'NumSolveOptions.Display = ''off'';\n');
%     fprintf(fid,'NumSolveOptions.MaxIterations = %f;\n',...
%             obj.NumSolveMaxIterations);
%     fprintf(fid,'NumSolveOptions.FunctionTolerance = %f;\n',...
%             obj.NumSolvePrecision);
%     fprintf(fid,'NumSolveOptions.OptimalityTolerance = %f;\n',...
%             obj.NumSolvePrecision);
%     fprintf(fid,'NumSolveOptions.StepTolerance = %f;\n',...
%             obj.NumSolvePrecision);
    fprintf(fid,['[NumSolveSolution,NumSolveResidual,NumSolveRC,',...
                     'NumSolveOutput] = ',...
                     'fsolve(@NumSolveEq,NumSolveGuess,NumSolveOptions);\n']);
    fprintf(fid,'Mats.NumSolveParamRC = NumSolveRC;\n');
%     fprintf(fid,'if NumSolveRC~=0\n');
    fprintf(fid,'if NumSolveRC~=1\n');
    fprintf(fid,'    Mats.Status = 0;\n');
    txt = 'NumSolveParam solution not normal.';
    fprintf(fid,'    Mats.StatusMessage = [Mats.StatusMessage,''%s''];\n',...
            txt);
    fprintf(fid,'    if op.verbose\n');
    fprintf(fid,'        fprintf(fid,''Warning: %s\\n'');\n',txt);
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'%s = NumSolveSolution(%.0f);\n',...
                obj.NumSolveParam.Names{j},j);
    end
end

if obj.CompoundParam.N>0
    fprintf(fid,'\n%% Map compound parameters\n');
    if obj.NumSolveParam.N>0
        fprintf(fid,'function EvalCompoundParam \n');
        txt = '    ';
    else
        txt = '';
    end
    for j=1:obj.CompoundParam.N
        fprintf(fid,'%s%s = %s;\n',txt,obj.CompoundParam.Names{j},...
                obj.CompoundExpressions{j});
    end
    if obj.NumSolveParam.N>0
        fprintf(fid,'end \n');
    end
end

% Combine NumSolve, and Compound into AuxParam
obj.AuxParam.Names = [obj.NumSolveParam.Names;obj.CompoundParam.Names];
obj.AuxParam.PrettyNames = [obj.NumSolveParam.PrettyNames;
                    obj.CompoundParam.PrettyNames];

fprintf(fid,'if op.StoreParam\n');
fprintf(fid,'    Mats.AuxParam = nan(%.0f,1);\n',obj.AuxParam.N);
for j=1:obj.AuxParam.N
    fprintf(fid,'    Mats.AuxParam(%.0f) = %s;\n',j,obj.AuxParam.Names{j});
%     fprintf(fid,'    Mats.AuxParam.%1$s = %1$s;\n',...
%             obj.AuxParam.Names{j});
end
fprintf(fid,'end\n');

if obj.ObsVar.N>0
    fprintf(fid,'\n%% Observation equations\n');
    H0 = -jacobian(ObsEq,ObsVar_t);
    H = jacobian(ObsEq,StateVar_t);
    HBar = simplify(ObsEq-H*StateVar_t.'+H0*ObsVar_t.');
%     SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
    SymMats.ObsEq.HBar = H0\HBar;
    SymMats.ObsEq.H = H0\H;
    idxEq = ( any(jacobian(ObsEq,StateVar_tF)~=0,2) ...
            & any(jacobian(ObsEq,StateVar_tL)~=0,2) ...
            & any(jacobian(ObsEq,ShockVar_t)~=0,2) );
    if any(idxEq)
        fprintf(2,'Equations violating model structure rules:\n');
        fprintf(2,'Obs Eq #%.0f\n',find(idxEq));
        fclose(fid);
        error('Cannot have leads, lags, or shocks in Obs equations.')
    end
    MatNames = fieldnames(SymMats.ObsEq);
    nCols = [1,obj.StateVar.N];
    fprintf(fid,'if op.StoreObsEq || op.StoreKF\n');
%     obj.ObsEqMats = struct;
    for jM=1:length(MatNames)
%         Mj = MatNames{jM};
%         nMj = obj.ObsVar.N*nCols(jM);
%         fh = cell(nMj,1);
%         for j=1:nMj
%             fh{j} = matlabFunction(SymMats.ObsEq.(Mj)(j),...
%                                     'Vars',SymParamAll);
%         end
%         obj.ObsEqMats.(MatNames{jM}) = @(x)buildmat(fh,x,obj.ObsVar.N,nCols(jM));
        fprintf(fid,'    ObsEq.%s = [...\n',MatNames{jM});
        for jeq=1:obj.ObsVar.N
            fprintf(fid,'       ');
            for jc=1:nCols(jM)
                fprintf(fid,' %s',...
                        char(eval(sprintf('SymMats.ObsEq.%s(jeq,jc)',...
                                          MatNames{jM}))));
                if jc==nCols(jM)
                    fprintf(fid,';\n');
                else
                    fprintf(fid,',');
                end
            end
        end
        fprintf(fid,'        ];\n');
    end
    fprintf(fid,'end\n');
    fprintf(fid,'if op.StoreObsEq\n');
    fprintf(fid,'    Mats.ObsEq = ObsEq;\n');
    fprintf(fid,'end\n');
end

fprintf(fid,'\n%% State equation matrices\n');
SymMats.StateEq.GammaBar = [];
SymMats.StateEq.Gamma0 = -jacobian(StateEq,StateVar_tF);
SymMats.StateEq.Gamma1 = jacobian(StateEq,StateVar_t);
SymMats.StateEq.Gamma4 = jacobian(StateEq,StateVar_tL);
SymMats.StateEq.Gamma2 = jacobian(StateEq,ShockVar_t);
% SymMats.StateEq.GammaBar = jacobian(StateEq,one);
SymMats.StateEq.GammaBar = simplify(StateEq ...
    + SymMats.StateEq.Gamma0*StateVar_tF.' ...
    - SymMats.StateEq.Gamma1*StateVar_t.' ...
    - SymMats.StateEq.Gamma4*StateVar_tL.' ...
    - SymMats.StateEq.Gamma2*ShockVar_t.');
idxEq = ( any(SymMats.StateEq.Gamma0~=0,2) ...
          & any(SymMats.StateEq.Gamma4~=0,2) );
if any(idxEq)
    fprintf(2,'Equations violating model structure rules:\n');
    fprintf(2,'State Eq #%.0f\n',find(idxEq));
    fclose(fid);
    error('Cannot have both leads and lags in same State equation.')
end

MatNames = fieldnames(SymMats.StateEq);
nCols = [1,obj.StateVar.N,obj.StateVar.N,obj.StateVar.N,obj.ShockVar.N];
for jM=1:length(MatNames)
    fprintf(fid,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:obj.StateVar.N
        fprintf(fid,'   ');
        for jc=1:nCols(jM)
            fprintf(fid,' %s',...
                    char(eval(sprintf('SymMats.StateEq.%s(jeq,jc)',...
                                      MatNames{jM}))));
            if jc==nCols(jM)
                fprintf(fid,';\n');
            else
                fprintf(fid,',');
            end
        end
    end
    fprintf(fid,'    ];\n\n');
end
fprintf(fid,'StateEq.Gamma3 = eye(%.0f);\n\n',obj.StateVar.N);
fprintf(fid,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',...
        obj.StateVar.N);
fprintf(fid,'StateEq.Gamma0(cv,:) = -StateEq.Gamma1(cv,:);\n');
fprintf(fid,'StateEq.Gamma1(cv,:) = StateEq.Gamma4(cv,:);\n');
fprintf(fid,'StateEq.Gamma3(:,cv) = [];\n');
fprintf(fid,'if ~all(all(StateEq.Gamma4(~cv,:)==0,2))\n');
fprintf(fid,'    error(''Incorrect system reduction'')\n');
fprintf(fid,'end\n\n');
fprintf(fid,'StateEq = rmfield(StateEq,''Gamma4'');\n');
fprintf(fid,'if op.StoreStateEq\n');
fprintf(fid,'    Mats.StateEq = StateEq;\n');
fprintf(fid,'end\n');

fprintf(fid,'\n%% Solve REE\n');
fprintf(fid,'if op.SolveREE\n');
fprintf(fid,...
        '    [REE,fmat,fwt,ywt,gev] = solveree(StateEq,...\n');
fprintf(fid,...
        '        op.GensysAuthor,op.fid,op.verbose,op.gensys{:});\n');
fprintf(fid,'    Mats.REE = REE;\n');
fprintf(fid,'    if ~all(REE.eu==1);\n');
fprintf(fid,'        Mats.Status = 0;\n');
fprintf(fid,['        Mats.StatusMessage = [Mats.StatusMessage,''REE ' ...
                 'solution not normal.''];\n']);
fprintf(fid,'    end\n');
fprintf(fid,'end\n');

if obj.ObsVar.N>0
    fprintf(fid,'\n%% Kalman Filter matrices\n');
    fprintf(fid,'if op.SolveREE && op.StoreKF\n');
    fprintf(fid,'    if all(Mats.REE.GBar(:)==0)\n');
    fprintf(fid,'        KF.StateVarBar = zeros(%.0f,1);\n',...
            obj.StateVar.N);
    fprintf(fid,'    else\n');
    fprintf(fid,...
            '        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
            obj.StateVar.N);
    fprintf(fid,'    end\n');
    fprintf(fid,...
            '    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

    if ~isempty(obj.KFInitState)
        fprintf(fid,'    s00 = [...\n');
        for jeq=1:obj.StateVar.N
            fprintf(fid,'        %.16f;\n',obj.KFInitState(jeq));
        end
        fprintf(fid,'        ];\n\n');
    else
        fprintf(fid,'    KF.s00 = zeros(%.0f,1);\n\n',obj.StateVar.N);
    end

    if ~isempty(obj.KFInitVariance)
        fprintf(fid,'    sig00 = [...\n');
        for jeq=1:obj.StateVar.N
            fprintf(fid,'       ');
            for jc=1:obj.StateVar.N
                fprintf(fid,' %0.16f',obj.KFInitVariance(jeq,jc));
                if jc==obj.StateVar.N
                    fprintf(fid,';\n');
                else
                    fprintf(fid,',');
                end
            end
        end
        fprintf(fid,'    ];\n\n');
        fprintf(fid,'    sig00rc = 0;\n');
    else
        fprintf(fid,...
                '    [sig00,sig00rc] = lyapcsd(REE.G1,REE.G2*REE.G2'');\n');
        fprintf(fid,...
                '    sig00 = real(sig00); sig00 = (sig00+sig00'')/2;\n');
        fprintf(fid,'    if sig00rc~=0\n');
        fprintf(fid,'        Mats.Status = 0;\n');
        txt = 'Could not find unconditional variance.';
        fprintf(fid,...
                '        Mats.StatusMessage = [Mats.StatusMessage,''%s''];\n',...
                txt);
        fprintf(fid,'        if op.verbose\n');
        fprintf(fid,...
                '            fprintf(fid,''Warning: %s\\n'');\n',txt);
        fprintf(fid,'        end\n\n');
        fprintf(fid,'    end\n\n');
    end
    fprintf(fid,'    KF.sig00 = sig00;\n');
    fprintf(fid,'    KF.sig00rc = sig00rc;\n');
    fprintf(fid,'    Mats.KF = KF;\n');
    fprintf(fid,'end\n');
end

if obj.AuxVar.N>0
    fprintf(fid,'\n%% Auxiliary equations matrices\n');
    SymMats.AuxEq.PhiBar = [];
    SymMats.AuxEq.Phi1 = jacobian(AuxEq,StateVar_t);
    SymMats.AuxEq.Phi2 = jacobian(AuxEq,ShockVar_t);
    SymMats.AuxEq.Phi3 = jacobian(AuxEq,StateVar_tF);
    SymMats.AuxEq.Phi4 = jacobian(AuxEq,StateVar_tL);
    SymMats.AuxEq.PhiBar = simplify(AuxEq ...
        - SymMats.AuxEq.Phi1*StateVar_t.' ...
        - SymMats.AuxEq.Phi2*ShockVar_t.' ...
        - SymMats.AuxEq.Phi3*StateVar_tF.' ...
        - SymMats.AuxEq.Phi4*StateVar_tL.');
    MatNames = {'PhiBar','Phi1','Phi2','Phi3','Phi4'};
    nCols = [1,obj.StateVar.N,obj.ShockVar.N,obj.StateVar.N,obj.StateVar.N];
    fprintf(fid,'if op.StoreAuxEq || op.StoreAuxREE\n');
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        fprintf(fid,'    AuxEq.%s = [...\n',Mj);
        for jeq=1:obj.AuxVar.N
            fprintf(fid,'       ');
            for jc=1:nCols(jM)
                fprintf(fid,' %s',...
                        char(eval(sprintf('SymMats.AuxEq.%s(jeq,jc)',Mj))));
                if jc==nCols(jM)
                    fprintf(fid,';\n');
                else
                    fprintf(fid,',');
                end
            end
        end
        fprintf(fid,'        ];\n\n');
    end
    fprintf(fid,'    if op.StoreAuxEq\n');
    fprintf(fid,'        Mats.AuxEq = AuxEq;\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'    if op.SolveREE && op.StoreAuxREE\n');
    fprintf(fid,'        if ~isempty(REE.G1)\n');
    fprintf(fid,['            AuxREE.GBar = ',...
                     'AuxEq.PhiBar+AuxEq.Phi3*REE.GBar',...
                     '+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.GBar',...
                     ';\n']);
    fprintf(fid,['            AuxREE.G1 = AuxEq.Phi4',...
                     '+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G1;\n']);
    fprintf(fid,['            AuxREE.G2 = AuxEq.Phi2',...
                     '+(AuxEq.Phi1+AuxEq.Phi3*REE.G1)*REE.G2;\n']);
    fprintf(fid,'        else\n');
    fprintf(fid,'            AuxREE.GBar = [];\n');
    fprintf(fid,'            AuxREE.G1 = [];\n');
    fprintf(fid,'            AuxREE.G2 = [];\n');
    fprintf(fid,'        end\n');
    fprintf(fid,'        Mats.AuxREE = AuxREE;\n');
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
end

% close file
fprintf(fid,'end\n');
fclose(fid);

%% Save handle to function
obj.mats = str2func(MatsFN);

%% Test function
Mats = obj.mats(obj.Param.Values);
if Mats.Status==0
    fprintf('Warning: REE solution not normal for Param.Values.\n')
end

%% Save timer
obj.TimeTracker.stop(ttName)
