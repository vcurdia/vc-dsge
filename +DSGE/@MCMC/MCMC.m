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
        nchains = 4;
        x0 = [];
        jumpscale = 2.4;
        jumpvar = [];
        jumpnewvarweight = 1;
        keeplogs = 1;
    end % properties
    
    methods
        function obj = MCMC(model)
            if nargin>0
                obj.model = model;
                obj.jumpvar = model.Post.Var;
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



