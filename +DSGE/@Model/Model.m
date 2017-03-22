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
        Param = initiatenames;
        NumSolveParam = initiatenames;
        NumSolveEq
        CompoundParam = initiatenames;
        AuxParam = initiatenames;
        ObsVar = initiatenames;
        StateVar = initiatenames;
        ShockVar = initiatenames;
        AuxVar = initiatenames;
        ObsEq
        StateEq
        AuxEq
        KFInitState
        KFInitVariance
        GensysAuthor = 'CS';
        NumSolvePrecision = 1e-6;
        NumSolveMaxIterations = 500;
        TimeElapsed
    end
    
    properties (SetAccess = protected)
        mats
    end
   
    methods
        function obj = Model(Name)
            if nargin>0
                obj.TimeElapsed = TimeTracker;
                obj.Name = Name;
            end
        end
        
        function set.Param(obj,p)
            if isstruct(p)
                obj.Param = p;
            else
                obj.Param = setnames('Param',p);
                if obj.Param.N>0
                    obj.Param.Values = [p{:,2}]';
                end
            end
        end
        
        function set.NumSolveParam(obj,p)
            if isstruct(p)
                obj.NumSolveParam = p;
            else
                obj.NumSolveParam = setnames('NumSolveParam',p);
                if obj.NumSolveParam.N>0
                    obj.NumSolveParam.Guess = [p{:,2}]';
                end
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
            if isstruct(p)
                obj.CompoundParam = p;
            else
                obj.CompoundParam = setnames('CompoundParam',p);
                if obj.CompoundParam.N>0
                    obj.CompoundParam.Expressions = p(:,2);
                end
            end
        end
    
        function set.ObsVar(obj,p)
            if isstruct(p)
                obj.ObsVar = p;
            else
                obj.ObsVar = setnames('ObsVar',p);
            end
        end
    
        function set.StateVar(obj,p)
            if isstruct(p)
                obj.StateVar = p;
            else
                obj.StateVar = setnames('StateVar',p);
            end
        end
    
        function set.ShockVar(obj,p)
            if isstruct(p)
                obj.ShockVar = p;
            else
                obj.ShockVar = setnames('ShockVar',p);
            end
        end
    
        function set.AuxVar(obj,p)
            if isstruct(p)
                obj.AuxVar = p;
            else
                obj.AuxVar = setnames('AuxVar',p);
                if obj.AuxVar.N>0
                    obj.AuxEq = p(:,2);
                end
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
        
    end %methods
    
end %class

function pp = setnames(pType,p)
    if ~ismember(pType,{'Param','NumSolveParam','CompoundParam',...
                        'ObsVar','StateVar','ShockVar','AuxVar'})
        error('SetProperty called with invalid type.')
    end
    pp = struct;
    [pp.N,nc] = size(p);
    pp.Names = p(:,1);
    if nc==(2+... 
            ismember(pType,{'Param','NumSolveParam','CompoundParam','AuxVar'}))
        pp.PrettyNames = p(:,nc);
    else
        pp.PrettyNames = pp.Names;
    end
end
        
function checkeq(obj,eqType)
    nEq = length(obj.([eqType,'Eq']));
    nVar = obj.([eqType,'Var']).N;
    if nEq~=nVar
        error(['Number of %sEq (%i) does not match number of ' ...
               'variables (%i).'],eqType,nEq,nVar)
    end
end

function s = initiatenames
    s.Names = {};
    s.PrettyNames = {};
    s.N = 0;
end


