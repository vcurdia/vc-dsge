classdef Param < matlab.mixin.Copyable

% DSGE.Param class
% 
% DSGE object representing parameters, used in the Model object.
%
% See also:
% DSGE.Model
%
% Created: November 7, 2016
% Copyright 2016-2018 Vasco Curdia
    
    properties 
% Names List of parameter names to be used in equations and commands.
        Names 

% PrettyNames - List of parameter names formatted for LaTeX
%  Use this property to create LaTeX output and figures.
        PrettyNames

% Values Vector of parameter values used in calibrated simulations.
        Values
        
% PriorDist String representing prior distribution type.
        PriorDist
        
% PriorMean Prior mean.
        PriorMean
        
% PriorSD Prior standard deviation.
        PriorSD
    end
    
    properties (SetAccess = protected)
        %N Number of parameters in instance
        N = 0;
    end
    
    methods
        
        function obj = Param(p)
            if nargin>0
                [np,nc] = size(p);
                obj.Names = p(:,1);
                if ismember(nc,[2,3])
                    obj.Values = [p{:,2}]';
                elseif nc>3
                    obj.PriorDist = p(:,2);
                    obj.PriorMean = [p{:,3}]';
                    for j=1:np
                        if isempty(p{j,4})
                            p{j,4} = 0;
                        end
                    end
                    obj.PriorSD = [p{:,4}]';
                    obj.Values = obj.PriorMean;
                end
                if ismember(nc,[3,5])
                    obj.PrettyNames = p(:,nc);
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
        
        function showvalues(obj,x)
            if nargin<2
                x = obj.Values;
            end
            v = struct;
            for j=1:obj.N
                v.(obj.Names{j}) = x(j);
            end
            disp(v)
        end
    
    end %methods
    
end %class
