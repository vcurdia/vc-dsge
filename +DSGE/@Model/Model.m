classdef Model
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
        FileName = struct;
        PlotDir = struct;
        TimeElapsed = struct;
        Param = InitiateNames;
        FixParam = InitiateNames;
        NumSolveParam = InitiateNames;
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
    end
   
    methods
        function obj = Model(Name)
            if nargin>0
                obj.Name = Name;
            end
        end
        
        function obj = set.Param(obj,p)
            if isstruct(p)
                obj.Param = p;
            else
                obj.Param = SetNames('Param',p);
                if obj.Param.N>0
                    obj.Param.PriorDist = p(:,2);
                    obj.Param.PriorMean = [p{:,3}]';
                    obj.Param.PriorSD = [p{:,4}]';
                end
            end
        end
        
        function obj = set.FixParam(obj,p)
            if isstruct(p)
                obj.FixParam = p;
            else
                obj.FixParam = SetNames('FixParam',p);
                if obj.FixParam.N>0
                    obj.FixParam.Values = [p{:,2}]';
                end
            end
        end
    
        function obj = set.NumSolveParam(obj,p)
            if isstruct(p)
                obj.NumSolveParam = p;
            else
                if length(p)<2
                    error(['Need to specify NumSolveParam as cell array ',...
                           'with two elements.'])
                end
                eq = p{2};
                p = p{1};
                obj.NumSolveParam = SetNames('NumSolveParam',p);
                if obj.NumSolveParam.N>0
                    obj.NumSolveParam.Guess = [p{:,2}]';
                    if isempty(eq) || (length(eq)<obj.NumSolveParam.N)
                        error(['Not enough equations specified for ',...
                               'NumSolveParam.'])
                    else
                        obj.NumSolveParam.Eq = eq;
                    end
                end
            end
        end
    
        function obj = set.CompoundParam(obj,p)
            if isstruct(p)
                obj.CompoundParam = p;
            else
                obj.CompoundParam = SetNames('CompoundParam',p);
                if obj.CompoundParam.N>0
                    obj.CompoundParam.Expressions = p(:,2);
                end
            end
        end
    
        function obj = set.ObsVar(obj,p)
            if isstruct(p)
                obj.ObsVar = p;
            else
                obj.ObsVar = SetNames('ObsVar',p);
            end
        end
    
        function obj = set.StateVar(obj,p)
            if isstruct(p)
                obj.StateVar = p;
            else
                obj.StateVar = SetNames('StateVar',p);
            end
        end
    
        function obj = set.ShockVar(obj,p)
            if isstruct(p)
                obj.ShockVar = p;
            else
                obj.ShockVar = SetNames('ShockVar',p);
            end
        end
    
        function obj = set.AuxVar(obj,p)
            if isstruct(p)
                obj.AuxVar = p;
            else
                obj.AuxVar = SetNames('AuxVar',p);
                if obj.AuxVar.N>0
                    obj.AuxEq = p(:,2);
                end
            end
        end
        
        function obj = set.ObsEq(obj,eq)
            obj.ObsEq = eq;
            CheckEq(obj,'Obs')
        end
    
        function obj = set.StateEq(obj,eq)
            obj.StateEq = eq;
            CheckEq(obj,'State')
        end
        
        function out = Mats(obj,x,varargin)
            out = feval(obj.FileName.Mats,[x;obj.FixParam.Values],varargin{:});
        end
        
    end %methods
    
end %class

function pp = SetNames(pType,p)
    if ~ismember(pType,{...
        'Param','FixParam','NumSolveParam','CompoundParam',...
        'ObsVar','StateVar','ShockVar','AuxVar',...
                       })
        error('SetProperty called with invalid type.')
    end
    pp = struct;
    [pp.N,nc] = size(p);
    pp.Names = p(:,1);
    if nc==(2+... 
            ismember(pType,{'FixParam','NumSolveParam','CompoundParam',...
                            'AuxVar'})+...
            3*strcmp(pType,'Param'))
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


