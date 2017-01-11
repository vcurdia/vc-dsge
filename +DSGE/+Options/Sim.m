classdef Sim < handle
% DSGE.Options.Sim class
% 
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: January 8, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 by Vasco Curdia
    
    properties
        Model
        FigPanels
        Test
    end
   
    methods
        
        function obj = Sim(m,varargin)
            obj.Model = m;
            obj = UpdateOptions(obj,varargin{:});
        end
        
    end %methods
    
end %class

