classdef Model < matlab.mixin.Copyable

% DSGE.Model class
%
% This is the main object of the DSGE. It describes all the variables,
% equations, and parameters of the DSGE model andincludes methods to solve and
% simulate the DSGE.
%
% Variables
% ---------
% 
%   ObsVar (optional) 
%   Observation variables. Their names need to match with variable names in the 
%   DSGE.Data object.
%
%   StateVar
%   State variables, including jump variables and pre-determined variables, 
%   endogenous or exogenous. It can also include any additional auxiliary 
%   variables needed to satisfy equation canonic form, such as adding lags or 
%   leads beyond one period.
% 
%   ShockVar
%   Innovations to the exogenous variables. Assumed to be iid normal 
%   distributed.
%
%   AuxVar (optional) 
%   Variables not needed to solve REE but which are useful to track for later 
%   simulations. Any variable that is not strictly needed for REE solution 
%   should be included in this object, so that the state space is kept as 
%   small as possible to improve performance.
%
% Variable objects do not include any time subscripts of any sort, just the
% variable names. In the equations need to reference variables always with a
% time subscript, subject to the following conventions:
%   x_t refers to x(t)
%   x_tF refers to E_t[x(t+1)] 
%   x_tL refers to x(t-1)
% 
% Equations
% ---------
%
%   ObsEq (optional)
%   Observation equations. Link ObsVar to StateVar:
%     0 = HBar + H*StateVar_t - H0*ObsVar_t
%   rules:
%     - no leads or lags fo any variables
%     - no ShockVar_t or AuxVar_t
%     
%   StateEq
%   State equations with laws of motion of economy in Chris Sims gensys 
%   canonical form:
%     0 = GammaBar + Gamma1*StateVar_tL + Gamma2*ShockVar_t + Gamma3*eta_t
%         - Gamma0*StateVar_t
%   where eta_t is an endogenous expectation error. The codes identify 
%   equations with forward looking components and automatically rearrange
%   matrices to fit in this canonical form.
%   rules:
%     - cannot have both leads and lags in same equation
%     - in order to use leads and lags in same equation create artificial 
%       variables, e.g. xL_t = x_tL means that new variable 'xL' is the lag 
%       of 'x'
%     - for higher order leads or lags use auxiliary variables as needed.
%   The solution of the REE yields:
%     StateVar_t = REE.GBar + REE.G1*StateVar_tL + REE.G2*ShockVar_t
%   
%   AuxEq (optional) 
%   Link AuxVar to StateVar:
%     AuxVar_t = PhiBar + Phi1*StateVar_t + Phi2*ShockVar_t 
%                + Phi3*StateVar_tF + Phi4*StateVar_tL
%   After solving REE and plugging into this equation we get
%     AuxVar_t = AuxREE.GBar + AuxREE.G1*StateVar_tL + AuxREE.G2*ShockVar_t
%
% Parameters
% ----------
% 
% See also:
% setupMyDSGE, DSGE.Var, DSGE.Param, solveree, gensysvb
%
% Created: November 7, 2016
% Copyright 2016-2017 Vasco Curdia
    
    properties
        Name = '';
        Param = DSGE.Param;
        NumSolveParam = DSGE.Param;
        NumSolveEq
        CompoundParam = DSGE.Param;
        CompoundExpressions
        AuxParam = DSGE.Param;
        ObsVar = DSGE.Var;
        StateVar = DSGE.Var;
        ShockVar = DSGE.Var;
        AuxVar = DSGE.Var;
        ObsEq
        StateEq
        AuxEq
        KFInitState
        KFInitVariance
        GensysAuthor = 'CS';
        NumSolvePrecision = 1e-6;
        NumSolveMaxIterations = 500;
        TimeTracker
        mats
    end
    
    methods
        function obj = Model(name)
            obj.TimeTracker = TimeTracker;
            if nargin>0
                fprintf('\n*** Preparing model\n')
                fprintf('%s\n',name)
                obj.Name = name;
            end
        end
        
        function set.Param(obj,p)
            if iscell(p)
                obj.Param = DSGE.Param(p);
            else
                obj.Param = p;
            end
        end
        
        function set.NumSolveParam(obj,p)
            if iscell(p)
                obj.NumSolveParam = DSGE.Param(p);
            else
                obj.NumSolveParam = p;
            end
        end
        
        function set.NumSolveEq(obj,eq)
            obj.NumSolveEq = eq;
            neq = length(eq);
            if neq~=obj.NumSolveParam.N
                error(['Number of NumSolveEq (%i) does not match number of ',...
                       'NumSolveParam (%i).'],neq,obj.NumSolveParam.N)
            end
        end
    
        function set.CompoundParam(obj,p)
            if iscell(p)
                obj.CompoundExpressions = p(:,2);
                p(:,2) = {nan};
                obj.CompoundParam = DSGE.Param(p);
            else
                obj.CompoundParam = p;
            end
        end
    
        function set.ObsVar(obj,v)
            if iscell(v)
                obj.ObsVar = DSGE.Var(v);
            else
                obj.ObsVar = v;
            end
        end
    
        function set.StateVar(obj,v)
            if iscell(v)
                obj.StateVar = DSGE.Var(v);
            else
                obj.StateVar = v;
            end
        end
    
        function set.ShockVar(obj,v)
            if iscell(v)
                obj.ShockVar = DSGE.Var(v);
            else
                obj.ShockVar = v;
            end
        end
    
        function set.AuxVar(obj,v)
            if iscell(v)
                obj.AuxEq = v(:,2);
                v(:,2) = [];
                obj.AuxVar = DSGE.Var(v);
            else
                obj.AuxVar = v;
            end
        end
        
        function set.ObsEq(obj,eq)
            obj.ObsEq = eq;
            checkeq(obj,'Obs')
        end
    
        function set.StateEq(obj,eq)
            obj.StateEq = eq;
            checkeq(obj,'State')
        end
        
        function showparamvalues(obj)
            p = struct;
            for j=1:obj.Param.N
                p.(obj.Param.Names{j}) = obj.Param.Values(j);
            end
            fprintf('Parameter values:\n')
            disp(p)
        end
        
        function setparamvalues(obj,ParNames,ParValues)
            if ischar(ParNames), ParNames = {ParNames}; end
            np = length(ParNames);
            nv = length(ParValues);
            if np>0 && np==nv
                [tf,idxp] = ismember(ParNames,obj.Param.Names);
                obj.Param.Values(idxp(tf)) = ParValues(tf);
                if ~all(tf)
                    fprintf('Invalid parameter names ignored:\n')
                    fprintf('  %s\n',ParNames{~tf})
                end
            else
                error(['List of parameter names is empty or has different ' ...
                       'length from list of values.'])
            end
        end
        
%         function Mats = evalmats(obj,matname,x)
%             Mats = struct;
%             matfields = fieldnames(obj.(matname));
%             for j=1:length(matfields)
%                 Mats.(matfields{j}) = obj.(matname).(matfields{j})(x);
%             end
%         end
            
    end %methods
    
    methods(Static)
    
    end %staticmethods
    
end %class

function checkeq(obj,eqType)
    nEq = length(obj.([eqType,'Eq']));
    nVar = obj.([eqType,'Var']).N;
    if nEq~=nVar
        error(['Number of %sEq (%i) does not match number of ' ...
               'variables (%i).'],eqType,nEq,nVar)
    end
end


