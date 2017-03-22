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
        ParamNames
        Dist
        Mean
        SD
        TimeElapsed = TimeTracker;
    end
   
    properties (SetAccess = protected)
        NParam
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
    end
    
    methods
        function obj = Prior(model,p)
            if nargin>1 && ~isempty(p)
                [obj.NParam,nc] = size(p);
                model.Param.N = obj.NParam;
                obj.ParamNames = p(:,1);
                model.Param.Names = obj.ParamNames;
                if nc==5
                    model.Param.PrettyNames = p(:,nc);
                else
                    model.Param.PrettyNames = model.Param.Names;
                end
                obj.Dist = p(:,2);
                obj.Mean = [p{:,3}]';
                for j=1:obj.NParam
                    if isempty(p{j,4})
                        p{j,4} = 0;
                    end
                end
                obj.SD = [p{:,4}]';
                model.Param.Values = obj.Mean;
                model.Param.EstimateIdx = ~ismember(obj.Dist,{'C'});
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

