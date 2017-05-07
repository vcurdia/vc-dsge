classdef Var

% DSGE.Var class
% 
% DSGE object representing variables, used in the model object.
% 
% See also:
% setupMyDSGE, DSGE.Model
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
% In not explicitly specified, then Names is used.
        PrettyNames 
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
        
        function obj = set.Names(obj,names)
            obj.Names = names;
            obj.N = length(names);
            if length(obj.PrettyNames)~=obj.N
                obj.PrettyNames = names;
            end
        end
        
        function obj = set.PrettyNames(obj,prettynames)
            if length(prettynames)==obj.N
                obj.PrettyNames = prettynames;
            else
                error('Length of PrettyNames must match number of variables.')
            end
        end
        
    end %methods
    
end %class
