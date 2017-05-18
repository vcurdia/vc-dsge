classdef Prior < matlab.mixin.Copyable
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
        Median
        Prc05
        Prc95
        Sample
        DistParam
        PDFCmd
        RndCmd
        LPDFCorrection
        TimeTracker
    end
   
    methods
        function obj = Prior(m,varargin)
            if nargin>0
                fprintf('\n*** Preparing prior\n')
                obj.TimeTracker = TimeTracker;
                obj.Model = m;
                if m.Param.N==0
                    error('There are no parameters in the model.')
                end
                obj.Dist = m.Param.PriorDist;
                obj.Mean = m.Param.PriorMean;
                obj.SD = m.Param.PriorSD;
                obj.analyzedist
                obj.analyzeparam(varargin{:})
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
        
    end %methods
    
end %class

