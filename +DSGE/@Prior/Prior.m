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
        PDFCmd
        RndCmd
        LPDFCorrection
        TimeElapsed
    end
    
    methods
        function obj = Prior(m,p)
            if nargin>0
                fprintf('\n*** Initiating prior\n')
                obj.TimeElapsed = TimeTracker;
                obj.Model = m;
            end
            if nargin>1 && ~isempty(p)
                [np,nc] = size(p);
                if nc>3
                    obj.Dist = p(:,2);
                    obj.Mean = [p{:,3}]';
                    for j=1:np
                        if isempty(p{j,4})
                            p{j,4} = 0;
                        end
                    end
                    obj.SD = [p{:,4}]';
                    p(:,[2,4]) = [];
                else
                    obj.Dist(1:np,1) = {'C'};
                    obj.SD(1:np,1) = 0;
                end
                m.Param = p;
                obj.EstimateIdx = ~ismember(obj.Dist,{'C'});
                obj.analyzedist
            end
        end
        
        function set.Model(obj,m)
            if isempty(obj.Dist) || ( m.Param.N==length(obj.Dist) )
                obj.Model = m;
            else
                error(['Cannot link prior to model with different number of ' ...
                       'parameters.'])
            end
        end
        
        function xd = draw(obj,nDraws)
            if nargin<2 || isempty(nDraws)
                nDraws = 1; 
            end
            xd = nan(obj.Model.Param.N,nDraws);
            for j=1:obj.Model.Param.N
                xd(j,:) = obj.RndCmd{j}(nDraws);
            end
        end
        
        function p = pdf(obj,x)
            p = repmat(exp(obj.LPDFCorrection),1,size(x,2));
            for j=1:obj.Model.Param.N
                p = p.*obj.PDFCmd{j}(x(j,:));
            end
        end
        
        function p = lpdf(obj,x)
            p = log(obj.pdf(x));
        end
        
        function new = copy(obj)
            new = DSGE.Prior;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for j = 1:length(pList)
                new.(pList{j}) = obj.(pList{j});
            end
        end
    
    end %methods
    
end %class

