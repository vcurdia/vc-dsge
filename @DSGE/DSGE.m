classdef DSGE
% DSGE class
% 
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: November 7, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia
    
    properties
        Spec = '';
        SpecPath = './';
        FileName = struct;
%         PlotDir = struct;
%        Options = struct;
        TimeElapsed = struct;
        Param
        FixParam
        NumSolveParam
        AuxParam
        ObsVar
        StateVar
        ShockVar
        AuxVar
        ObsEq
        StateEq
        AuxEq
        NParam = 0;
        NFixParam = 0;
        NNumSolveParam = 0;
        NAuxParam = 0;
        NObsVar = 0;
        NStateVar = 0;
        NShockVar = 0;
        NAuxVar = 0;
        GensysAuthor = 'CS';
    end
   
    methods
        function obj = DSGE(Spec,SpecPath)
            if nargin>0
                obj.Spec = Spec;
            end
            if nargin<2
                if ~isempty(obj.Spec)
                    obj.SpecPath = [Spec,'/'];
                end
            end
        end
        
        function obj = set.SpecPath(obj,SpecPath)
            if ~ismember(SpecPath,{'','./'})
                obj.SpecPath = SpecPath;
                mkdir(obj.SpecPath)
            end
        end
        
        function obj = set.Param(obj,p)
            pType = 'Param';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
            if obj.NParam>0
                obj.Param.PriorDist = p(:,2);
                obj.Param.PriorMean = [p{:,3}]';
                obj.Param.PriorSD = [p{:,4}]';
            end
        end
        
        function obj = set.FixParam(obj,p)
            pType = 'FixParam';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
            if obj.NFixParam>0
                obj.FixParam.Values = [p{:,2}]';
            end
        end
    
        function obj = set.NumSolveParam(obj,p)
            if length(p)<2
                error(['Need to specify NumSolveParam as cell array with two ' ...
                       'elements.'])
            end
            eq = p{2};
            p = p{1};
            pType = 'NumSolveParam';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
            if obj.NNumSolveParam>0
                obj.NumSolveParam.Guess = [p{:,2}]';
                if isempty(eq) || (length(eq)<obj.NNumSolveParam)
                    error('Not enough equations specified for NumSolveParam.')
                else
                    obj.NumSolveParam.Eq = eq;
                end
            end
        end
    
        function obj = set.AuxParam(obj,p)
            pType = 'AuxParam';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
            if obj.NAuxParam>0
                obj.AuxParam.Expressions = p(:,2);
            end
        end
    
        function obj = set.ObsVar(obj,p)
            pType = 'ObsVar';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
        end
    
        function obj = set.StateVar(obj,p)
            pType = 'StateVar';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
        end
    
        function obj = set.ShockVar(obj,p)
            pType = 'ShockVar';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
        end
    
        function obj = set.AuxVar(obj,p)
            pType = 'AuxVar';
            [obj.(['N',pType]),obj.(pType)] = SetNames(pType,p);
            if obj.NAuxVar>0
                obj.AuxEq = p(:,2);
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
    
    end
    
end

function [np,pp] = SetNames(pType,p)
    if ~ismember(pType,{...
        'Param','FixParam','NumSolveParam','AuxParam',...
        'ObsVar','StateVar','ShockVar','AuxVar',...
                       })
        error('SetProperty called with invalid type.')
    end
    pp = struct;
    [np,nc] = size(p);
    if np==0, return, end
    pp.Names = p(:,1);
    if nc==(2+... 
            ismember(pType,{'FixParam','NumSolveParam','AuxParam','AuxVar'})+...
            3*strcmp(pType,'Param'))
        pp.PrettyNames = p(:,nc);
    else
        pp.PrettyNames = pp.Names;
    end
end
        
function CheckEq(obj,eqType)
    nEq = length(obj.([eqType,'Eq']));
    nVar = obj.(['N',eqType,'Var']);
    if nEq~=nVar
        error(['Number of %sEq (%i) does not match number of ' ...
               'variables (%i).'],eqType,nEq,nVar)
    end
end

