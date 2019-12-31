classdef Var < matlab.mixin.Copyable

% DSGE.Var class
% 
% DSGE object representing variables, used in the model object.
% 
% See also:
% setupDSGE, DSGE.Model
%
% Created: November 7, 2016 
% Copyright 2016-2017 Vasco Curdia
    
    properties
% Names List of variable names to be used in equations and commands
%
% These are the names that will show up in equations and model manipulations.
        Names 
        
% PrettyNames List of variable names formatted for LaTeX output
%
% These are expressions that will represent the variables when plotting or 
% creating LaTeX output.
% 
% If not explicitly specified, then Names is used.
        PrettyNames 
        
% Scale
% 
% scale for plotting. Default is 1 for every variable.
        Scale
    end
    
    properties (SetAccess = protected)
% N number of variables in instance
%
% This property is automatically populated every time that Names is changed.
        N = 0; 
    end
    
    methods
        
        function obj = Var(v)
            if nargin>0
                [nv,nc] = size(v);
                obj.Names = v(:,1);
                if nc>1
                    obj.PrettyNames = v(:,2);
                end
            end
        end
        
        function set.Names(obj,names)
            obj.Names = names;
            obj.N = length(names);
            if length(obj.PrettyNames)~=obj.N
                obj.PrettyNames = names;
            end
            obj.Scale = ones(obj.N,1);
        end
        
        function set.PrettyNames(obj,prettynames)
            if length(prettynames)==obj.N
                obj.PrettyNames = prettynames;
            else
                error('Length of PrettyNames must match number of variables.')
            end
        end
        
        function add(obj,v)
            if ~iscell(v)
                error('Variables to add need to be in cell array.')
            end
            [nv,nc] = size(v);
            names = [obj.Names;v(:,1)];
            prettynames = [obj.PrettyNames;v(:,nc)];
            scale = [obj.Scale;ones(nv,1)];
            obj.Names = names;
            obj.PrettyNames = prettynames;
            obj.Scale = scale;
        end

        function setscale(obj,v)
            [tf,idx] = ismember(v(:,1),obj.Names);
            obj.Scale(idx) = [v{:,2}];
        end

        function v1 = merge(obj,varargin)
            v1 = copy(obj);
            isDuplicates = false;
            for j=1:nargin-1
                v = varargin{j};
                if ~strcmp(class(v),'DSGE.Var')
                    error(['Cannot merge Var. Input needs to be instance of ' ...
                           'DSGE.Var'])
                end
                if v.N==0, continue, end
                if v1.N==0
                    v1 = copy(v);
                    continue
                end
                tf = ismember(v.Names,v1.Names);
                isDuplicates = (isDuplicates || any(tf));
                if ~all(tf)
                    names = [v1.Names;v.Names(~tf)];
                    prettynames = [v1.PrettyNames;v.PrettyNames(~tf)];
                    scale = [v1.Scale;v.Scale(~tf)];
                    v1.Names = names;
                    v1.PrettyNames = prettynames;
                    v1.Scale = scale;
                end
            end
%             if isDuplicates
%                 fprintf('Duplicate variable names found.\n')
%             end
        end
        
        function [v1,idx] = subset(obj,names)
            [tf,idx] = obj.ismember(names);
            idx = idx(tf);
            v1 = DSGE.Var([obj.Names(idx),obj.PrettyNames(idx)]);
            v1.Scale = obj.Scale(idx);
        end
        
        function [tf,idx] = ismember(obj,names)
            if isa(names,'DSGE.Var'), names = names.Names; end
            [tf,idx] = ismember(names,obj.Names);
        end
        
        function prettynames=findprettynames(obj,names)
            [tf,idx] = ismember(names,obj.Names);
            if ~all(tf)
                error('Could not find %s\n',names{~tf})
            end
            prettynames = obj.PrettyNames(idx);
        end
        
    end %methods
    
end %class
