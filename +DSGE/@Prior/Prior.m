classdef Prior < handle
% DSGE.Prior class
% 
% See also:
% SetupMyDSGE, DSGE.Model
%
% ...........................................................................
%
% Created: March 19, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 Vasco Curdia
    
    properties
        Model
        Dist
        Mean
        SD
        Mode
        Sample
    end
   
    properties (SetAccess = protected)
    end
    
    methods
        function obj = Prior(m,p)
            if nargin>0
                obj.Model = m;
                m.Prior = obj;
            end
            if nargin>1 && ~isempty(p)
                [np,nc] = size(p);
                m.Param.N = np;
                m.Param.Names = p(:,1);
                if nc==5
                    m.Param.PrettyNames = p(:,nc);
                else
                    m.Param.PrettyNames = m.Param.Names;
                end
                if ismember(nc,[2,3])
                    m.Param.Values = [p{:,2}]';
                else
                    obj.Dist = p(:,2);
                    obj.Mean = [p{:,3}]';
                    for j=1:np
                        if isempty(p{j,4})
                            p{j,4} = 0;
                        end
                    end
                    obj.SD = [p{:,4}]';
                    m.Param.Values = obj.Mean;
                    m.Param.EstimateIdx = ~ismember(obj.Dist,{'C'});
                end
            end
        end
        
    end %methods
    
end %class

