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
    end
   
    properties (SetAccess = protected)
        EstimateIdx
        Mode
        Median
        Prc05
        Prc95
        Sample
        DistParam
        LPdfCmd
        PdfCmd
        RndCmd
        LPdfCorrection
        TimeElapsed
    end
    
    methods
        function obj = Prior(model,p)
            if nargin>0
                obj.Model = model;
            end
            if nargin>1 && ~isempty(p)
                obj.TimeElapsed = TimeTracker;
                [np,nc] = size(p,2);
                model.Param.Names = p(:,1);
                if nc==5
                    model.Param.PrettyNames = p(:,nc);
                end
                obj.Dist = p(:,2);
                obj.Mean = [p{:,3}]';
                for j=1:np
                    if isempty(p{j,4})
                        p{j,4} = 0;
                    end
                end
                obj.SD = [p{:,4}]';
                model.Param.Values = obj.Mean;
                obj.EstimateIdx = ~ismember(obj.Dist,{'C'});
                obj.analyzedist
            end
        end
        
        function xd = draw(obj,nDraws)
            if nargin<2 || isempty(nDraws)
                nDraws = 1; 
            end
            xd = nan(obj.NParam,nDraws);
            for j=1:obj.NParam
                xd(j,:) = obj.RndCmd{j}(nDraws);
            end
        end
        
    end %methods
    
end %class

