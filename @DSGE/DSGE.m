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
        Name = '';
        Path = './';
        FileName = struct;
        Report
%         PlotDir = struct;
        TimeElapsed = struct;
        Param
        FixParam
        NumSolveParam
        CompoundParam
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
        NCompoundParam = 0;
        NAuxParam = 0;
        NObsVar = 0;
        NStateVar = 0;
        NShockVar = 0;
        NAuxVar = 0;
        GensysAuthor = 'CS';
        KFInit
        Prior
        Post
    end
   
    methods
        function obj = DSGE(Name,Path)
            if nargin>0
                obj.Name = Name;
            end
            if nargin<2
                if ~isempty(obj.Name)
                    obj.Path = [Name,'/'];
                end
            end
        end
        
        function obj = set.Path(obj,Path)
            obj.Path = Path;
            if ~ismember(Path,{'','./'}) && ~isdir(obj.Path)
                mkdir(obj.Path)
            end
        end
        
        function obj = set.Param(obj,p)
            if isstruct(p)
                obj.Param = p;
            else
                pType = 'Param';
                [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
                if obj.NParam>0
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
                pType = 'FixParam';
                [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
                if obj.NFixParam>0
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
                pType = 'NumSolveParam';
                [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
                if obj.NNumSolveParam>0
                    obj.NumSolveParam.Guess = [p{:,2}]';
                    if isempty(eq) || (length(eq)<obj.NNumSolveParam)
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
                pType = 'CompoundParam';
                [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
                if obj.NCompoundParam>0
                    obj.CompoundParam.Expressions = p(:,2);
                end
            end
        end
    
        function obj = set.ObsVar(obj,p)
            pType = 'ObsVar';
            [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
        end
    
        function obj = set.StateVar(obj,p)
            pType = 'StateVar';
            [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
        end
    
        function obj = set.ShockVar(obj,p)
            pType = 'ShockVar';
            [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
        end
    
        function obj = set.AuxVar(obj,p)
            pType = 'AuxVar';
            [obj.(['N',pType]),obj.(pType)] = setnames(pType,p);
            if obj.NAuxVar>0
                obj.AuxEq = p(:,2);
            end
        end
        
        function obj = set.ObsEq(obj,eq)
            obj.ObsEq = eq;
            checkeq(obj,'Obs')
        end
    
        function obj = set.StateEq(obj,eq)
            obj.StateEq = eq;
            checkeq(obj,'State')
        end
        
        function out = mats(obj,x,varargin)
            out = feval(obj.FileName.Mats,[x;obj.FixParam.Values],varargin{:});
        end
    
        function out = drawprior(obj,varargin)
            out = feval(obj.FileName.DrawPrior,varargin{:});
        end
    
        function obj = tracktime(obj,action,start)
            if start
                obj.TimeElapsed.(action) = tic;
            else
                obj.TimeElapsed.(action) = toc(obj.TimeElapsed.(action));
                fprintf('\n%s %s\n\n',action,vctoc(obj.TimeElapsed.(action)))
            end
        end
        
    end
    
end

function [np,pp] = setnames(pType,p)
    if ~ismember(pType,{...
        'Param','FixParam','NumSolveParam','CompoundParam',...
        'ObsVar','StateVar','ShockVar','AuxVar',...
                       })
        error('SetProperty called with invalid type.')
    end
    pp = struct;
    [np,nc] = size(p);
    if np==0, return, end
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
        
function checkeq(obj,eqType)
    nEq = length(obj.([eqType,'Eq']));
    nVar = obj.(['N',eqType,'Var']);
    if nEq~=nVar
        error(['Number of %sEq (%i) does not match number of ' ...
               'variables (%i).'],eqType,nEq,nVar)
    end
end

