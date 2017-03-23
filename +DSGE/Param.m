classdef Param
% DSGE.Param class
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
        Values
    end
    
    properties (SetAccess = protected)
        N = 0;
    end
    
    methods
        
        function obj = Param(p)
            if nargin>0
                [np,nc] = size(p);
                obj.Names = p(:,1);
                if nc>1
                    obj.Values = [p{:,2}]';
                end
                if nc>2
                    obj.PrettyNames = p(:,3);
                end
            end
        end
        
        function obj = set.Names(obj,names)
            obj.Names = names;
            obj.N = length(names);
            if length(obj.PrettyNames)~=obj.N
                obj.PrettyNames = names;
            end
            if length(obj.Values)~=obj.N
                obj.Values = nan(obj.N,1);
            end
        end
        
        function obj = set.PrettyNames(obj,prettynames)
            if length(prettynames)==obj.N
                obj.PrettyNames = prettynames;
            else
                error('Length of PrettyNames must match number of parameters.')
            end
        end
        
        function obj = set.Values(obj,values)
            if length(values)==obj.N
                obj.Values = values;
            else
                error('Length of Values must match number of parameters.')
            end
        end            
        
    end %methods
    
end %class
