classdef Var
% DSGE.Var class
% 
% See also:
% setupMyDSGE, DSGE.Model
%
% ...........................................................................
%
% Created: November 7, 2016 by Vasco Curdia
% 
% Copyright (C) 2016-2017 Vasco Curdia
    
    properties
        Names
        PrettyNames
    end
    
    properties (SetAccess = protected)
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
