function genmatsobj(obj,matspath)

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
fprintf('\n*** Generating DSGE mats\n')
ttName = 'GenMatsObj';
obj.TimeTracker.start(ttName)

%% checks
if (obj.StateVar.N==0) || isempty(obj.StateEq)
    error('Cannot proceed without specifying state variables and equations')
end

%% Sym Params
list = {'Param','NumSolveParam','CompoundParam'};
symlist = struct;
for j=1:length(list)
    jstr = list{j};
    if obj.(jstr).N>0, vcsym(obj.(jstr).Names{:}), end
    symlist.(jstr) = sym(zeros(1,obj.(jstr).N));
    for jp=1:obj.(jstr).N
        symlist.(jstr)(jp) = eval(obj.(jstr).Names{jp});
    end
end

% Combine NumSolve, and Compound into AuxParam
obj.AuxParam.Names = [obj.NumSolveParam.Names;obj.CompoundParam.Names];
obj.AuxParam.PrettyNames = [obj.NumSolveParam.PrettyNames;
                    obj.CompoundParam.PrettyNames];

% Define AllParam
AllParam = [obj.Param.Names;obj.AuxParam.Names];
symAllParam = [symlist.Param,symlist.NumSolveParam,symlist.CompoundParam];
nAllParam = length(AllParam);


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


%% function to evaluate compound parameters
ftmp = cell(obj.CompoundParam.N,1);
fj = sym(zeros(obj.CompoundParam.N,1));
for j=1:obj.CompoundParam.N
    fj(j) = eval(obj.CompoundExpressions{j});
    for jj=j-1:-1:1
        fj(j) = subs(fj(j),symlist.CompoundParam(jj),fj(jj));
    end
    ftmp{j} = matlabFunction(fj(j),'Vars',[symlist.Param,symlist.NumSolveParam]);
end
obj.MatFcn.CompoundParam = @(x)buildmat(ftmp,x,obj.CompoundParam.N,1);

%% function to evaluate NumSolveEq
ftmp = cell(obj.NumSolveParam.N,1);
for j=1:obj.NumSolveParam.N
    fj = eval(obj.NumSolveEq{j});
    ftmp{j} = matlabFunction(fj,'Vars',symAllParam);
end
obj.MatFcn.NumSolveEq = @(x)buildmat(ftmp,[x;obj.MatFcn.CompoundParam(x)],...
                      obj.NumSolveParam.N,1);

%% functions to evaluate ObsEq
if obj.ObsVar.N>0
    H0 = -jacobian(ObsEq,ObsVar_t);
    H = jacobian(ObsEq,StateVar_t);
    HBar = simplify(ObsEq-H*StateVar_t.'+H0*ObsVar_t.');
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
end

%% State equation matrices
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

%% Auxiliary equations matrices
if obj.AuxVar.N>0
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
end

listEq = {'ObsEq','StateEq','AuxEq'};
for jEq=1:length(listEq)
    Eqj = listEq{jEq};
    if strcmp(Eqj,'ObsEq')
        nRows = obj.ObsVar.N;
        nCols = [1,obj.StateVar.N];
    elseif strcmp(Eqj,'StateEq')
        nRows = obj.StateVar.N;
        nCols = [1,obj.StateVar.N,obj.StateVar.N,obj.StateVar.N, ...
                 obj.ShockVar.N];
    elseif strcmp(Eqj,'AuxEq')
        nRows = obj.AuxVar.N;
        nCols = [1,obj.StateVar.N,obj.ShockVar.N,obj.StateVar.N, ...
                 obj.StateVar.N];
    end
    if nRows==0,continue,end
    MatNames = fieldnames(SymMats.(Eqj));
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        nMj = nRows*nCols(jM);
        fj = cell(nMj,1);
        for j=1:nMj
            fj{j} = matlabFunction(SymMats.(Eqj).(Mj)(j),'Vars',symAllParam);
        end
        obj.MatFcn.(Eqj).(MatNames{jM}) = @(x)buildmat(fj,x,nRows,nCols(jM));
    end
end

%% Test function
mats = obj.evalmats(obj.Param.Values);
obj.AuxParam.Values = mats.AuxParam;
if mats.Status==0
    fprintf('Warning: REE solution not normal for Param.Values.\n')
end

%% Save timer
obj.TimeTracker.stop(ttName)
