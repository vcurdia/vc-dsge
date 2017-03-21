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
        TimeElapsed = TimeTracker;
    end
   
    properties (SetAccess = protected)
        Mode
        Median
        Prc05
        Prc95
        NParam
        Sample
        DistParam
        LPdfCmd
        PdfCmd
        RndCmd
        LPdfCorrection
    end
    
    methods
        function obj = Prior(m,p)
            if nargin>0
                obj.Model = m;
                m.Prior = obj;
            end
            if nargin>1 && ~isempty(p)
                [obj.NParam,nc] = size(p);
                m.Param.N = obj.NParam;
                m.Param.Names = p(:,1);
                if nc==5
                    m.Param.PrettyNames = p(:,nc);
                else
                    m.Param.PrettyNames = m.Param.Names;
                end
                obj.Dist = p(:,2);
                obj.Mean = [p{:,3}]';
                for j=1:obj.NParam
                    if isempty(p{j,4})
                        p{j,4} = 0;
                    end
                end
                obj.SD = [p{:,4}]';
                m.Param.Values = obj.Mean;
                m.Param.EstimateIdx = ~ismember(obj.Dist,{'C'});
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

