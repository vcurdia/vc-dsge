function genmats(obj,matspath)

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
% DSGE, setupMyDSGE
%
%
% Created: January 22, 2016 by Vasco Curdia
% Copyright (C) 2016-2017 Vasco Curdia


%% Preamble
fprintf('Generating DSGE mats\n')

%% checks
if nargin==1|| isempty(matspath)
    matspath = '';
else
    if ~strcmp(matspath(end),'/'), matspath = [matspath,'/']; end
    if ~isdir(matspath), mkdir(matspath), end
    addpath(matspath)
end

if (obj.StateVar.N==0) || isempty(obj.StateEq)
    error('Cannot proceed without specifying state variables and equations')
end

%% Sym Params
list = {'','NumSolve','Composite'};
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
% SymCompositeParam = sList.CompositeParam;
% SymParamAll = [SymParam,SymNumsolveParam,SymCompositeParam];

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
if obj.AuxVar.N>0
    nEq = length(obj.AuxEq);
    if nEq~=obj.AuxVar.N
        error('Number of AuxVar (%.0f) and AuxEq (%.0f) do not match.',...
              obj.AuxVar.N,nEq)
    end
    AuxEq = sym(zeros(obj.AuxVar.N,1));
    for j=1:obj.AuxVar.N
        vj = [obj.AuxVar.Names{j},'_t'];
        eval([vj,' = ',obj.AuxEq{j},';'])
        AuxEq(j) = eval(vj);
        eval([vj,'F = ',strrep(obj.AuxEq{j},'_t','_tF'),';'])
        eval([vj,'L = ',strrep(obj.AuxEq{j},'_t','_tL'),';'])
    end
end

%% Build Observation equations
if obj.ObsVar.N>0
    nEq = length(obj.ObsEq);
    if nEq~=obj.ObsVar.N
        error('Number of ObsVar (%.0f) and ObsEq (%.0f) do not match.',...
              obj.ObsVar.N,nEq)
    end
    ObsEq = sym(zeros(obj.ObsVar.N,1));
    for j=1:obj.ObsVar.N
        ObsEq(j) = eval(obj.ObsEq{j});
    end
end

%% Build State equations
nEq = length(obj.StateEq);
if nEq~=obj.StateVar.N
    error('Number of StateVar (%.0f) and StateEq (%.0f) do not match.',...
          obj.StateVar.N,nEq)
end
StateEq = sym(zeros(obj.StateVar.N,1));
for j=1:obj.StateVar.N
    StateEq(j) = eval(obj.StateEq{j});
end


%% Generate code to eval matrices

MatsFN = sprintf('%smats',obj.Name);

% Initiate file
fid = fopen([matspath,MatsFN,'.m'],'wt');
fprintf(fid,'function Mats = %s(x,varargin)\n\n',MatsFN);
fprintf(fid,'%% Created: %.0f/%.0f/%.0f %.0f:%.0f:%.0fs\n',clock);

fprintf(fid,'\n%% options\n');
fprintf(fid,'op.FID = 1;\n');
fprintf(fid,'op.Verbose = 0;\n');
fprintf(fid,'op = updateoptions(op,varargin{:});\n');

fprintf(fid,'\n%% Initiate Status\n');
fprintf(fid,'Mats.Status = 1;\n');
fprintf(fid,'Mats.StatusMessage = '''';\n');

fprintf(fid,'\n%% Map parameters\n');
for j=1:obj.Param.N
    fprintf(fid,'%s = x(%.0f);\n',obj.Param.Names{j},j);
end
fprintf(fid,'Mats.Param = x;\n');

if obj.CompositeParam.N>0 || obj.NumSolveParam.N>0
    fprintf(fid,'\n%% Initialize composite parameters\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'%s = [];\n',obj.NumSolveParam.Names{j});
    end
    for j=1:obj.CompositeParam.N
        fprintf(fid,'%s = [];\n',obj.CompositeParam.Names{j});
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
    fprintf(fid,'        EvalCompositeParam\n');
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
    fprintf(fid,'    if op.Verbose\n');
    fprintf(fid,'        fprintf(op.FID,''Warning: %s\\n'');\n',txt);
    fprintf(fid,'    end\n');
    fprintf(fid,'end\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fid,'%s = NumSolveSolution(%.0f);\n',...
                obj.NumSolveParam.Names{j},j);
    end
end

if obj.CompositeParam.N>0
    fprintf(fid,'\n%% Map composite parameters\n');
    if obj.NumSolveParam.N>0
        fprintf(fid,'function EvalCompositeParam \n');
        txt = '    ';
    else
        txt = '';
    end
    for j=1:obj.CompositeParam.N
        fprintf(fid,'%s%s = %s;\n',txt,obj.CompositeParam.Names{j},...
                obj.CompositeExpressions{j});
    end
    if obj.NumSolveParam.N>0
        fprintf(fid,'end \n');
    end
end

% Combine NumSolve, and Composite into AuxParam
obj.AuxParam.Names = [obj.NumSolveParam.Names;obj.CompositeParam.Names];
obj.AuxParam.PrettyNames = [obj.NumSolveParam.PrettyNames;
                    obj.CompositeParam.PrettyNames];

fprintf(fid,'Mats.AuxParam = nan(%.0f,1);\n',obj.AuxParam.N);
for j=1:obj.AuxParam.N
    fprintf(fid,'Mats.AuxParam(%.0f) = %s;\n',j,obj.AuxParam.Names{j});
end
fprintf(fid,'if ~all(isreal(Mats.AuxParam))\n');
fprintf(fid,'    Mats.Status = 0;\n');
fprintf(fid,'    Mats.StatusMessage = ''Auxiliary parameters with no-real values.'';\n');
fprintf(fid,'    return\n');
fprintf(fid,'end\n\n');

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
        fprintf(fid,'ObsEq.%s = [...\n',MatNames{jM});
        for jeq=1:obj.ObsVar.N
            fprintf(fid,'    ');
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
        fprintf(fid,'    ];\n');
    end
    fprintf(fid,'Mats.ObsEq = ObsEq;\n');
end

fprintf(fid,'\n%% State equation matrices\n');
SymMats.StateEq.GammaBar = [];
SymMats.StateEq.Gamma0 = -jacobian(StateEq,StateVar_tF);
SymMats.StateEq.Gamma1 = jacobian(StateEq,StateVar_t);
SymMats.StateEq.Gamma4 = jacobian(StateEq,StateVar_tL);
SymMats.StateEq.Gamma2 = jacobian(StateEq,ShockVar_t);
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
fprintf(fid,'Mats.StateEq = StateEq;\n');

if obj.AuxVar.N>0
    fprintf(fid,'\n%% Auxiliary equations matrices\n');
    SymMats.AuxEq.PhiBar = [];
    SymMats.AuxEq.Phi = jacobian(AuxEq,StateVar_t);
    SymMats.AuxEq.PhiBar = simplify(AuxEq - SymMats.AuxEq.Phi*StateVar_t.');
    MatNames = {'PhiBar','Phi'};
    nCols = [1,obj.StateVar.N,obj.ShockVar.N,obj.StateVar.N,obj.StateVar.N];
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        fprintf(fid,'AuxEq.%s = [...\n',Mj);
        for jeq=1:obj.AuxVar.N
            fprintf(fid,'    ');
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
        fprintf(fid,'    ];\n\n');
    end
    fprintf(fid,'Mats.AuxEq = AuxEq;\n');
end

% close file
fprintf(fid,'end\n');
fclose(fid);

%% Save handle to function
obj.mats = str2func(MatsFN);

%% Test function
mats = obj.solveree(obj.Param.Values,'Verbose',1);
obj.AuxParam.Values = mats.AuxParam;


end

