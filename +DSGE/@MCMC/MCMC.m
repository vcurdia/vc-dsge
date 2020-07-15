classdef MCMC < matlab.mixin.Copyable

% DSGE.MCMC class
%
% See also:
% ../Example/estimatedsge, DSGE.Model
%
% Created: July 15, 2020
% Copyright 2020 Vasco Curdia
    
    properties
        model = [];
        ndraws = 100000;
        ndrawskeep = 100000;
    end % properties
    
    methods
        function obj = MCMC(model)
            if nargin>0
                obj.model = model;
            end
        end
        
        function set.ndraws(obj,n)
            obj.ndraws = n;
            obj.ndrawskeep = min(n,obj.ndrawskeep);
        end
        
        function set.ndrawskeep(obj,n)
            obj.ndrawskeep = min(n,obj.ndraws);
        end
        
    end %methods
    
    methods(Static)
    
    end %staticmethods
    
end %class



