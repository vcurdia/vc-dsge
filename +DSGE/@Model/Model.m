classdef Model < handle
% DSGE.Model class
% 
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: November 7, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia
    
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
            
        function new = copy(obj)
            new = DSGE.Model;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for j = 1:length(pList)
                new.(pList{j}) = obj.(pList{j});
            end
        end
    
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


