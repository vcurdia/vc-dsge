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
        Param = InitiateNames;
        NumSolveParam = InitiateNames;
        NumSolveEq
        CompoundParam = InitiateNames;
        AuxParam = InitiateNames;
        ObsVar = InitiateNames;
        StateVar = InitiateNames;
        ShockVar = InitiateNames;
        AuxVar = InitiateNames;
        ObsEq
        StateEq
        AuxEq
        KFInitState
        KFInitVariance
        GensysAuthor = 'CS';
        NumSolvePrecision = 1e-6;
        NumSolveMaxIterations = 500;
        Mats
        MatsFN
    end
   
    methods
        function obj = Model(Name)
            if nargin>0
                obj.Name = Name;
            end
        end
        
        function set.Param(obj,p)
            if isstruct(p)
                obj.Param = p;
            else
                obj.Param = SetNames('Param',p);
                if obj.Param.N>0
                    obj.Param.Values = [p{:,2}]';
                end
            end
        end
        
        function ShowParam(obj)
            p = struct;
            for j=1:obj.Param.N
                p.(obj.Param.Names{j}) = obj.Param.Values(j);
            end
            disp(p)
        end
        
        function set.NumSolveParam(obj,p)
            if isstruct(p)
                obj.NumSolveParam = p;
            else
                obj.NumSolveParam = SetNames('NumSolveParam',p);
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
                obj.CompoundParam = SetNames('CompoundParam',p);
                if obj.CompoundParam.N>0
                    obj.CompoundParam.Expressions = p(:,2);
                end
            end
        end
    
        function set.ObsVar(obj,p)
            if isstruct(p)
                obj.ObsVar = p;
            else
                obj.ObsVar = SetNames('ObsVar',p);
            end
        end
    
        function set.StateVar(obj,p)
            if isstruct(p)
                obj.StateVar = p;
            else
                obj.StateVar = SetNames('StateVar',p);
            end
        end
    
        function set.ShockVar(obj,p)
            if isstruct(p)
                obj.ShockVar = p;
            else
                obj.ShockVar = SetNames('ShockVar',p);
            end
        end
    
        function set.AuxVar(obj,p)
            if isstruct(p)
                obj.AuxVar = p;
            else
                obj.AuxVar = SetNames('AuxVar',p);
                if obj.AuxVar.N>0
                    obj.AuxEq = p(:,2);
                end
            end
        end
        
        function set.ObsEq(obj,eq)
            obj.ObsEq = eq;
            CheckEq(obj,'Obs')
        end
    
        function set.StateEq(obj,eq)
            obj.StateEq = eq;
            CheckEq(obj,'State')
        end
        
    end %methods
    
end %class

function pp = SetNames(pType,p)
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
        
function CheckEq(obj,eqType)
    nEq = length(obj.([eqType,'Eq']));
    nVar = obj.([eqType,'Var']).N;
    if nEq~=nVar
        error(['Number of %sEq (%i) does not match number of ' ...
               'variables (%i).'],eqType,nEq,nVar)
    end
end

function s = InitiateNames
    s.N = 0;
    s.Names = {};
    s.PrettyNames = {};
end


