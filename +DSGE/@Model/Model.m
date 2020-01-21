classdef Model < matlab.mixin.Copyable

% DSGE.Model class
%
% This is the main object of the DSGE. It describes all the variables,
% equations, and parameters of the DSGE model and includes methods to solve and
% simulate the DSGE.
%
% Click on the links for additional information on each property or method.
%
% Model Properties:
%
% * Properties describing model variables
%   ObsVar   - Observation variables (optional)
%   StateVar - State variables
%   ShockVar - Shock variables
%   AuxVar   - Auxiliary Variables (optional) 
% 
%   Variable objects do not include any time subscripts of any sort, just the
%   variable names. In the equations need to reference variables always with a
%   time subscript, subject to the following conventions:
%     x_t  refers to x(t)
%     x_tF refers to E_t[x(t+1)] 
%     x_tL refers to x(t-1)
% 
%
% * Properties describing model equations
%   ObsEq   - Observation equations (optional)
%   StateEq - State equations
%   AuxEq   - Auxiliary equations (optional) 
%
%   Assumed Model structure:
%
%     0 = HBar + H*StateVar_t - H0*ObsVar_t
%
%     0 = GammaBar + Gamma1*StateVar_t + Gamma2*ShockVar_t + Gamma4*StateVar_tL
%         - Gamma0*StateVar_tF
%
%     AuxVar_t = PhiBar + Phi*StateVar_t
%  
%   Note: cannot include both leads and lags in the same equation. If the model 
%         specification has equations with both leads and lags, then need to 
%         create artificial variables as needed. The example illustrates how to
%         accomplish that.
%
%   Method genmats will prepare all matrices and convert the state equations 
%   into the gensys' canonical form:
%     0 = GammaBar + Gamma1*StateVar_tL + Gamma2*ShockVar_t + Gamma3*eta_t
%         - Gamma0*StateVar_t
%
%   Model solution:
%     StateVar_t = REE.GBar     + REE.G1*StateVar_tL   + REE.G2*ShockVar_t
%     ObsVar_t   = ObsEq.HBar   + ObsEq.H*StateVar_t
%     AuxVar_t   = AuxEq.PhiBar + AuxEq.Phi*StateVar_t
%
%
% * Properties describing model parameters
%   Param          - Set of parameters to be calibrated or estimated
%   NumSolveParam  - Parameters to be solved numberically (optional)
%   CompositeParam - Parameters that are combinations of Param (optional)
%   AuxParam       - combined NumSolveParam and CompositeParam (not defined by 
%                    user)
%
%   To define Param can proceed in one of two ways, depending on how model will
%   be used.
%
%   For calibrated model simulations:
%   Set Param with a cell array in which each row has the following:
%   - name of parameter
%   - value of parameter
%   - LaTeX representation of parameter (optional)
%
%   For estimating model:
%   Set Param with a cell array in which each row has the following:
%   - name of parameter
%   - prior distribution code within the following list
%     'C'   - Calibrated parameter
%     'N'   - Normal
%     'TN'  - Truncated Normal (truncated at zero)
%     'B'   - Beta 
%     'G'   - Gamma
%     'IG1' - Inverse Gamma type 1 (for standard deviation parameters)
%     'IG2' - Inverse Gamma type 2 (for variances parameters)
%   - prior mean
%   - prior SD
%   - LaTeX representation of parameter (optional)
%   
% Refer to setupdsge in the example folder for a concrete case on how to setup
% the model, prior and posterior.
%
% 
% See also:
% ../Example/setupdsge, DSGE.Var, DSGE.Param, solveree, gensysvb
%
% Created: November 7, 2016
% Copyright 2016-2018 Vasco Curdia
    
    properties
% Name of the model specification
        Name = '';
        
% Param - Set of parameters to be calibrated or estimated
%
%   Main set of parameters set or estimated for a model specification.
% 
%   Can be set directly for model calibration or through DSGE.Prior to set it 
%   up for model estimation.
        Param = DSGE.Param;

% NumSolveParam - Parameters to be solved numberically (optional)
%
%   Parameters that need to be solved numerically, as functions of Param and 
%   CompositeParam.
%
%   System of nonlinear equations are set in model property NumSolveEq.
        NumSolveParam = DSGE.Param;

% NumSolveEq - system of equations to solve for NumSolveParam
%
%   Array with expressions for the system of nonlinear equations to be solved 
%   for NunmSolveParam.
        NumSolveEq
        
% CompositeParam - Parameters that are combinations of Param (optional)
%
%   Parameters that are combinations of the main Param. Used to renormalize 
%   parameters, de4fined ratios of parameters to use in equation expressions, 
%   or simply reference combinations of parameters to keep track (e.g. as 
%   checks to model behavior)
%
%   Their expressions are defined in model propery CompositeExpressions.
        CompositeParam = DSGE.Param;
        
% CompositeExpressions - Definitions of CompositeParam
%
%   Array of expressions with definitions of CompositeParam. 
        CompositeExpressions

% AuxParam - combined NumSolveParam and CompositeParam (not defined by user)
%
%   Combines all parameters set in NumSolveParam and CompounParam to be used in 
%   simulations. This way there is only one set of auxiliary parameters to keep 
%   track, rather than two.
        AuxParam = DSGE.Param;
        
% ObsVar - Observation variables (optional)
%   DSGE.Var object representing observation variables.
%   Their names need to match with variable names in the DSGE.Data object.
        ObsVar = DSGE.Var;
        
% StateVar - State variables
%   Include jump variables and pre-determined variables, endogenous or 
%   exogenous. It can also include any additional auxiliary variables needed to 
%   satisfy equation canonic form, such as adding lags or leads beyond one 
%   period.
        StateVar = DSGE.Var;
        
%  ShockVar - Shock variables
%   Innovations to the exogenous variables. Assumed to be iid normal 
%   distributed.
        ShockVar = DSGE.Var;
        
% AuxVar - Auxiliary Variables (optional) 
%   Variables not needed to solve REE but which are useful to track for later 
%   simulations. Restriction is that only depends on current period states and
%   a constant.
        AuxVar = DSGE.Var;

% ObsEq - Observation equations (optional)
%
%   Array with expressions for the equations linking ObsVar to StateVar:
%     0 = HBar + H*StateVar_t - H0*ObsVar_t
%   rules:
%     - no leads or lags for any variables
%     - no ShockVar_t or AuxVar_t
        ObsEq

% StateEq - State equations
%
%   Array with expressions for the laws of motion of economy in Chris Sims 
%   gensys canonical form:
%     0 = GammaBar + Gamma1*StateVar_tL + Gamma2*ShockVar_t + Gamma3*eta_t
%         - Gamma0*StateVar_t
%   where eta_t is an endogenous expectation error. The codes identify 
%   equations with forward looking components and automatically rearrange
%   matrices to fit in this canonical form.
%
%   Rules:
%     - cannot have both leads and lags in same equation
%     - in order to use leads and lags in same equation create artificial 
%       variables, e.g. xL_t = x_tL means that new variable 'xL' is the lag 
%       of 'x'
%     - for higher order leads or lags use auxiliary variables as needed.
%
%   The solution of the REE yields:
%     StateVar_t = REE.GBar + REE.G1*StateVar_tL + REE.G2*ShockVar_t
        StateEq
        
% AuxEq - Auxiliary equations (optional) 
%
%   Array with expressions of equations defining AuxVar as functions of 
%   StateVar and ShockVar.
%
%   Assumed structure:
%     AuxVar_t = PhiBar + Phi*StateVar_t
%   rules:
%     - no leads or lags for any variables
%     - no ShockVar_t
        AuxEq
        
% structure containing Prior information, if a prior exists
        Prior
        
% Data for the model
        Data
        
% structure containing Posterior information, if a posterior exists
        Post
        
% Other model properties and options
        KFInitState
        KFInitVariance
        GensysOptions
%         NumSolvePrecision = 1e-6;
%         NumSolveMaxIterations = 500;

% handle to the mats evaluation function
        mats
        
    end
    
    methods
        function obj = Model(name)
            if nargin>0
                fprintf('Model: %s\n',name)
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
    
        function set.CompositeParam(obj,p)
            if iscell(p)
                obj.CompositeExpressions = p(:,2);
                p(:,2) = {nan};
                obj.CompositeParam = DSGE.Param(p);
            else
                obj.CompositeParam = p;
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
        
%         function set.ObsEq(obj,eq)
%             obj.ObsEq = eq;
%             checkeq(obj,'Obs')
%         end
    
%         function set.StateEq(obj,eq)
%             obj.StateEq = eq;
%             checkeq(obj,'State')
%         end

        function mats = setparamvalues(obj,p)
            x = obj.Param.getvalues(p);
            obj.Param.Values = x;
            mats = obj.solveree(obj.Param.Values);
            obj.AuxParam.Values = mats.AuxParam;
        end

        function checkeq(obj,eqType)
            nEq = length(obj.([eqType,'Eq']));
            nVar = obj.([eqType,'Var']).N;
            if nEq~=nVar
                fprintf(['Warning: Number of %sEq (%i) does not match number ' ...
                         'of variables (%i).\n'],eqType,nEq,nVar)
            end
        end
        
        function prettynames = findvarprettynames(obj,names)
            vv = merge(obj.ObsVar,obj.StateVar,obj.AuxVar);
            [tf,idx] = ismember(names,vv.Names);
            if ~all(tf)
                error('Could not find %s\n',names{~tf})
            end
            prettynames = vv.PrettyNames(idx);
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



