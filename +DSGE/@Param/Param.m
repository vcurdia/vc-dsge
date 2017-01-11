classdef Param < handle
% DSGE.Param class
% 
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: January 10, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 Vasco Curdia
    
    properties
        Model
        N = 0;
        Names = {};
        PrettyNames = {};
        Values
        PriorDist
        PriorMean
        PriorSD
    end
   
    methods
        function obj = Param(m,p)
            if nargin>0
                obj.Model = m;
                m.Param = obj;
            end
            if nargin>1 && ~isempty(p)
                [obj.N,nc] = size(p);
                obj.Names = p(:,1);
                if nc==5
                    obj.PrettyNames = p(:,nc);
                else
                    obj.PrettyNames = obj.Names;
                end
                if ismember(nc,[2,3])
                    obj.Values = [p{:,2}]';
                else
                    obj.PriorDist = p(:,2);
                    obj.PriorMean = p(:,3);
                    obj.PriorSD = p(:,4);
                    obj.Values = obj.PriorMean;
                end
                
            end
        end
        
        
    end %methods
    
end %class

