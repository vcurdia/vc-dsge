classdef Posterior < handle
% DSGE.Posterior class
% 
% See also:
% setupMyDSGE, DSGE.Model, DSGE.Prior
%
% ...........................................................................
%
% Created: March 22, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 Vasco Curdia
    
    properties
        Model
        Prior
        Data
        Mode
        Mean
        SD
        Median
        Prc05
        Prc95
        Sample
    end
   
    properties (SetAccess = protected)
        TimeElapsed
    end
    
    methods
        function obj = Posterior(model,prior,data)
            if nargin>0
                obj.TimeElapsed = TimeTracker;
                obj.Model = model;
                obj.Prior = prior;
                obj.Data = data;
            end
        end
        
        function xd = draw(obj,nDraws)
            if nargin<2 || isempty(nDraws)
                nDraws = 1; 
            end
            xd = nan(obj.NParam,nDraws);
            error('Posterior draw method is not ready yet.')
        end
        
    end %methods
    
end %class

