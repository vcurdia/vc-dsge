function obj = GenMats(obj)

% GenMats
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
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

action = 'GenMats';
obj = obj.TrackTime(action,1);

fprintf('\n*** Generate DSGE mats\n')

%% -------------------------------------------------------------------

%% Prepare parameters, variables and equations

fprintf('Generating symbolic variables and systems of equations\n')

%% basic check
if (obj.StateVar.N==0) || isempty(obj.StateEq)
    error('Cannot without specifying state variables and states')
end

%% Sym Params
list = {'','Fix','NumSolve','Compound'};
for j=1:length(list)
    jstr = [list{j},'Param'];
    if obj.(jstr).N>0, vcSym(obj.(jstr).Names{:}), end
end
if obj.FixParam.N==0
    obj.FixParam.Values = [];
end


%% Constant
vcSym('one')

%% Obs Var
ObsVar_t = sym(zeros(1,obj.ObsVar.N)); 
for j=1:obj.ObsVar.N
    vj = [obj.ObsVar.Names{j},'_t'];
    vcSym(vj)
    ObsVar_t(j) = eval(vj);
end

%% State Var
StateVar_t = sym(zeros(1,obj.StateVar.N)); 
StateVar_tF = sym(zeros(1,obj.StateVar.N)); 
StateVar_tL = sym(zeros(1,obj.StateVar.N)); 
for j=1:obj.StateVar.N
    vj = [obj.StateVar.Names{j},'_t'];
    vcSym(vj,[vj,'F'],[vj,'L'])
    StateVar_t(j) = eval(vj);
    StateVar_tF(j) = eval([vj,'F']);
    StateVar_tL(j) = eval([vj,'L']);
end

%% Shock Var
ShockVar_t = sym(zeros(1,obj.ShockVar.N)); 
for j=1:obj.ShockVar.N
    vj = [obj.ShockVar.Names{j},'_t'];
    vcSym(vj)
    ShockVar_t(j) = eval(vj);
end
obj.Shocks2Show = obj.ShockVar.Names;

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

%% Set default FigPanels
PanelList = {};
if obj.ObsVar.N>0, PanelList{end+1} = 'ObsVar';end
PanelList{end+1} = 'StateVar';
if obj.AuxVar.N>0, PanelList{end+1} = 'AuxVar';end
jP = 0;
for jList=1:length(PanelList)
    Listj = PanelList{jList};
    nPanelj = ceil(obj.(Listj).N/obj.FigPanelMaxVar);
    for j=1:nPanelj
        jP = jP+1;
        if obj.FigPanelMaxVar==1
            obj.FigPanels(jP).Var = obj.(Listj).Names(j);
            obj.FigPanels(jP).PrettyNames = obj.(Listj).PrettyNames(j);
            obj.FigPanels(jP).Title = sprintf('%s_%s',Listj,...
                                          obj.FigPanels(jP).Var{1});
        else
            if nPanelj>1
                obj.FigPanels(jP).Title = sprintf('%s_%.0f',Listj,j);
            else
                obj.FigPanels(jP).Title = sprintf('%s',Listj);
            end
            obj.FigPanels(jP).Var = obj.(Listj).Names(...
                (j-1)*obj.FigPanelMaxVar+1:min(j*obj.FigPanelMaxVar,...
                                           obj.(Listj).N));
            obj.FigPanels(jP).PrettyNames = obj.(Listj).PrettyNames(...
                (j-1)*obj.FigPanelMaxVar+1:min(j*obj.FigPanelMaxVar,...
                                           obj.(Listj).N));
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

obj.FileName.Mats = sprintf('%s_Mats',obj.Name);

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

fprintf(fidMats,'\n%% Map parameters\n');
for j=1:obj.Param.N
    fprintf(fidMats,'%s = x(%.0f);\n',obj.Param.Names{j},j);
end
for j=1:obj.FixParam.N
    fprintf(fidMats,'%s = x(%.0f);\n',obj.FixParam.Names{j},obj.Param.N+j);
end
fprintf(fidMats,'if op.StoreParam\n');
for j=1:obj.Param.N
    fprintf(fidMats,'    Mats.Param.%1$s = %1$s;\n',obj.Param.Names{j});
end
for j=1:obj.FixParam.N
    fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',...
            obj.FixParam.Names{j});
end
fprintf(fidMats,'end\n');

if obj.CompoundParam.N>0 || obj.NumSolveParam.N>0
    fprintf(fidMats,'\n%% Initialize compound parameters\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'%s = [];\n',obj.NumSolveParam.Names{j});
    end
    for j=1:obj.CompoundParam.N
        fprintf(fidMats,'%s = [];\n',obj.CompoundParam.Names{j});
    end
end

if obj.NumSolveParam.N>0
    fprintf(fidMats,'\n%% NumSolve parameters\n');
    fprintf(fidMats,'function f=NumSolveEq(x)\n');
    fprintf(fidMats,'    for jx=1:size(x,2)\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'        %s = x(%.0f,jx);\n',...
                obj.NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'        EvalCompoundParam\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'        f(%.0f,jx) = %s;\n',j,...
                obj.NumSolveParam.Eq{j});
    end
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    fprintf(fidMats,'NumSolveGuess = [...\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'    %.16f;\n',obj.NumSolveParam.Guess(j));
    end
    fprintf(fidMats,'    ];\n');
%     fprintf(fidMats,['[NumSolveSolution,NumSolveRC] = csolvevb(@NumSolveEq,' ...
%                      'NumSolveGuess,[],%e,%.0f);\n'],...
%             obj.NumPrecision,obj.NumSolveMaxIterations);
    fprintf(fidMats,['NumSolveOptions = optimoptions(@fsolve);\n']);
    fprintf(fidMats,'NumSolveOptions.Display = ''off'';\n');
%     fprintf(fidMats,'NumSolveOptions.MaxIterations = %f;\n',...
%             obj.NumSolveMaxIterations);
%     fprintf(fidMats,'NumSolveOptions.FunctionTolerance = %f;\n',...
%             obj.NumSolvePrecision);
%     fprintf(fidMats,'NumSolveOptions.OptimalityTolerance = %f;\n',...
%             obj.NumSolvePrecision);
%     fprintf(fidMats,'NumSolveOptions.StepTolerance = %f;\n',...
%             obj.NumSolvePrecision);
    fprintf(fidMats,['[NumSolveSolution,NumSolveResidual,NumSolveRC,',...
                     'NumSolveOutput] = ',...
                     'fsolve(@NumSolveEq,NumSolveGuess,NumSolveOptions);\n']);
    fprintf(fidMats,'Mats.NumSolveParamRC = NumSolveRC;\n');
%     fprintf(fidMats,'if NumSolveRC~=0\n');
    fprintf(fidMats,'if NumSolveRC~=1\n');
    fprintf(fidMats,'    Mats.Status = 0;\n');
    txt = 'NumSolveParam solution not normal.';
    fprintf(fidMats,'    Mats.StatusMessage = [Mats.StatusMessage,''%s''];\n',...
            txt);
    fprintf(fidMats,'    if op.verbose\n');
    fprintf(fidMats,'        fprintf(fid,''Warning: %s\\n'');\n',txt);
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,'end\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'%s = NumSolveSolution(%.0f);\n',...
                obj.NumSolveParam.Names{j},j);
    end
    fprintf(fidMats,'if op.StoreParam\n');
    for j=1:obj.NumSolveParam.N
        fprintf(fidMats,'    Mats.NumSolveParam.%1$s = %1$s;\n',...
                obj.NumSolveParam.Names{j});
    end
    fprintf(fidMats,'end\n');
end

if obj.CompoundParam.N>0
    fprintf(fidMats,'\n%% Map compound parameters\n');
    if obj.NumSolveParam.N>0
        fprintf(fidMats,'function EvalCompoundParam \n');
        txt = '    ';
    else
        txt = '';
    end
    for j=1:obj.CompoundParam.N
        fprintf(fidMats,'%s%s = %s;\n',txt,obj.CompoundParam.Names{j},...
                obj.CompoundParam.Expressions{j});
    end
    if obj.NumSolveParam.N>0
        fprintf(fidMats,'end \n');
    end
end

% Combine Fix, NumSolve, and Compound into AuxParam
obj.AuxParam.N = obj.FixParam.N + obj.NumSolveParam.N + obj.CompoundParam.N;
obj.AuxParam.Names = [obj.FixParam.Names;
                    obj.NumSolveParam.Names;obj.CompoundParam.Names];
obj.AuxParam.PrettyNames = [obj.FixParam.PrettyNames;
                    obj.NumSolveParam.PrettyNames;
                    obj.CompoundParam.PrettyNames];

fprintf(fidMats,'if op.StoreParam\n');
for j=1:obj.NumSolveParam.N
    fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',...
            obj.NumSolveParam.Names{j});
end
for j=1:obj.CompoundParam.N
    fprintf(fidMats,'    Mats.AuxParam.%1$s = %1$s;\n',...
            obj.CompoundParam.Names{j});
end
fprintf(fidMats,'end\n');

if obj.ObsVar.N>0
    fprintf(fidMats,'\n%% Observation equations\n');
    H0 = -jacobian(ObsEq,ObsVar_t);
    SymMats.ObsEq.HBar = H0\jacobian(ObsEq,one);
    SymMats.ObsEq.H = H0\jacobian(ObsEq,StateVar_t);
    MatNames = fieldnames(SymMats.ObsEq);
    nCols = [1,obj.StateVar.N];
    fprintf(fidMats,'if op.StoreObsEq || op.StoreKF\n');
    for jM=1:length(MatNames)
        fprintf(fidMats,'    ObsEq.%s = [...\n',MatNames{jM});
        for jeq=1:obj.ObsVar.N
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
idxFL = (any(SymMats.StateEq.Gamma0~=0,2) & any(SymMats.StateEq.Gamma4~=0,2));
if any(idxFL)
    fprintf(2,['The following equations have both leads (''_tF'') and ',...
               'lags (''_tL'')\n']);
    fprintf(2,'State Eq #%.0f\n',find(idxFL));
    fclose(fidMats);
    error('Cannot have both leads and lags in same equation.')
end

MatNames = fieldnames(SymMats.StateEq);
nCols = [1,obj.StateVar.N,obj.StateVar.N,obj.StateVar.N,obj.ShockVar.N];
for jM=1:length(MatNames)
    fprintf(fidMats,'StateEq.%s = [...\n',MatNames{jM});
    for jeq=1:obj.StateVar.N
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
fprintf(fidMats,'StateEq.Gamma3 = eye(%.0f);\n\n',obj.StateVar.N);
fprintf(fidMats,'cv = (all(StateEq.Gamma0(1:%.0f,:)==0,2)~=0);\n',...
        obj.StateVar.N);
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
        obj.GensysAuthor);
fprintf(fidMats,'    Mats.REE = REE;\n');
fprintf(fidMats,'    if ~all(REE.eu==1);\n');
fprintf(fidMats,'        Mats.Status = 0;\n');
fprintf(fidMats,['        Mats.StatusMessage = [Mats.StatusMessage,''REE ' ...
                 'solution not normal.''];\n']);
fprintf(fidMats,'    end\n');
fprintf(fidMats,'end\n');

if obj.ObsVar.N>0
    fprintf(fidMats,'\n%% Kalman Filter matrices\n');
    fprintf(fidMats,'if op.StoreKF\n');
    fprintf(fidMats,'    if all(Mats.REE.GBar(:)==0)\n');
    fprintf(fidMats,'        KF.StateVarBar = zeros(%.0f,1);\n',...
            obj.StateVar.N);
    fprintf(fidMats,'    else\n');
    fprintf(fidMats,...
            '        KF.StateVarBar = (eye(%.0f)-REE.G1)\\REE.GBar;\n',...
            obj.StateVar.N);
    fprintf(fidMats,'    end\n');
    fprintf(fidMats,...
            '    KF.ObsVarBar = ObsEq.HBar + ObsEq.H*KF.StateVarBar;\n\n');

    if ~isempty(obj.KFInitState)
        fprintf(fidMats,'    s00 = [...\n');
        for jeq=1:obj.StateVar.N
            fprintf(fidMats,'        %.16f;\n',obj.KFInitState(jeq));
        end
        fprintf(fidMats,'        ];\n\n');
    else
        fprintf(fidMats,'    KF.s00 = zeros(%.0f,1);\n\n',obj.StateVar.N);
    end

    if ~isempty(obj.KFInitVariance)
        fprintf(fidMats,'    sig00 = [...\n');
        for jeq=1:obj.StateVar.N
            fprintf(fidMats,'       ');
            for jc=1:obj.StateVar.N
                fprintf(fidMats,' %0.16f',obj.KFInitVariance(jeq,jc));
                if jc==obj.StateVar.N
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
        fprintf(fidMats,'        Mats.Status = 0;\n');
        txt = 'Could not find unconditional variance.';
        fprintf(fidMats,...
                '        Mats.StatusMessage = [Mats.StatusMessage,''%s''];\n',...
                txt);
        fprintf(fidMats,'        if op.verbose\n');
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

if obj.AuxVar.N>0
    fprintf(fidMats,'\n%% Auxiliary equations matrices\n');
    MatNames = {'one','StateVar_t','StateVar_tF','StateVar_tL','ShockVar_t'};
    nCols = [1,obj.StateVar.N,obj.StateVar.N,obj.StateVar.N,obj.ShockVar.N];
    fprintf(fidMats,'if op.StoreAuxEq || op.StoreAuxREE\n');
    for jM=1:length(MatNames)
        Mj = MatNames{jM};
        SymMats.AuxEq.(Mj) = jacobian(AuxEq,eval(Mj));
        fprintf(fidMats,'    AuxEq.%s = [...\n',Mj);
        for jeq=1:obj.AuxVar.N
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
obj = obj.TrackTime(action,0);


%% -------------------------------------------------------------------
