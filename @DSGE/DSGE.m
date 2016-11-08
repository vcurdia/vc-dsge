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
        Param = struct;
        FixParam = struct;
        NumSolveParam = struct;
        AuxParam = struct;
%         Var = struct;
%         Eq = struct;
        NParam = 0;
        NFixParam = 0;
        NNumSolveParam = 0;
        NAuxParam = 0;
    end
   
    methods
        function obj = DSGE(Spec,SpecPath)
            if nargin>0
                obj.Spec = Spec;
            else
                error('Need to specify DSGE Spec.')
            end
            if nargin<2
                obj.SpecPath = [Spec,'/'];
            end
            mkdir(obj.SpecPath)
        end
        
        function obj = set.Param(obj,p)
            [obj.NParam,obj.Param] = SetParam(obj,'Param',p);
            if obj.NParam>0
                obj.Param.PriorDist = p(:,2);
                obj.Param.PriorMean = [p{:,3}]';
                obj.Param.PriorSD = [p{:,4}]';
            end
        end
        
        function obj = set.FixParam(obj,p)
            [obj.NFixParam,obj.FixParam] = SetParam(obj,'FixParam',p);
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
            [obj.NNumSolveParam,obj.NumSolveParam] = ...
                SetParam(obj,'NumSolveParam',p);
            if obj.NNumSolveParam>0
                obj.NumSolveParam.Guess = [p{:,2}]';
                if isempty(eq) || (length(eq)<obj.NNumSolveParam)
                    error('Not enough equations specified for NumSolveParam.')
                else
                    obj.NumSolveParam.Eq = eq;
                end
            end
        end
    
    end
    
end

function [np,pp] = SetParam(obj,pType,p)
    if ~ismember(pType,...
                 {'Param','FixParam','NumSolveParam','AuxParam'})
        error('SetParam called with invalid parameter type.')
    end
    pp = struct;
    [np,nc] = size(p);
    if np==0, return, end
    pp.Names = p(:,1);
    if nc==3+2*strcmp(pType,'Param')
        pp.PrettyNames = p(:,nc);
    else
        pp.PrettyNames = p.(pType).Names;
    end
end
        
