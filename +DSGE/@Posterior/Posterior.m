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
        EstimateIdx
        NEstimate
        Mode
        LPDFMode
        Mean
        SD
        Var
        Corr
        Median
        Prc05
        Prc95
        LogMgLikelihood
        MCMCStage
        MCMCSample
    end
   
    properties (SetAccess = protected)
        TimeElapsed
    end
    
    methods
        function obj = Posterior(model,prior,data)
            if nargin>0
                fprintf('\n*** Preparing posterior\n')
                obj.TimeElapsed = TimeTracker;
                obj.Model = model;
                obj.Prior = prior;
                obj.Data = data;
                obj.EstimateIdx = ~ismember(prior.Dist,{'C'});
                obj.NEstimate = sum(obj.EstimateIdx);
                obj.Mode = model.Param.Values;
                obj.LPDFMode = obj.lpdf(obj.Mode);
                fprintf('Posterior log-pdf using Param.Values is %0.4f.\n',...
                        obj.LPDFMode);
                obj.Mean = prior.Mean;
                obj.SD = prior.SD;
                obj.Var = diag(prior.SD.^2);
                obj.Median = prior.Median;
                obj.Prc05 = prior.Prc05;
                obj.Prc95 = prior.Prc95;
            end
        end
        
        function set.Model(obj,m)
            if isempty(obj.Mean) || ( length(obj.Mean)==m.Param.N )
                obj.Model = m;
            else
                error(['Cannot link posterior to model with different number ' ...
                       'of parameters.'])
            end
        end
        
        function set.Prior(obj,p)
            if isempty(obj.Mode) || ( length(obj.Mean)==length(p.Mean) )
                obj.Prior = p;
            else
                error(['Cannot link posterior to prior with different number ' ...
                       'of parameters.'])
            end
        end
        
        function xx = expandparam(obj,x)
            xx = repmat(obj.Model.Param.Values,1,size(x,2),size(x,3));
            xx(obj.EstimateIdx,:,:) = x;
        end
        
        function xd = draw(obj,nDraws)
            if nargin<2 || isempty(nDraws)
                nDraws = 1; 
            end
            if isempty(obj.MCMCSample.FileNameRedux)
                obj.mcmcredux
            end
            draws = load(obj.MCMCSample.FileNameRedux);
            xd = obj.expandparam(...
                draws.Param(:,randi(draws.N,1,nDraws)));
        end
        
        function new = copy(obj)
            new = DSGE.Posterior;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for j = 1:length(pList)
                new.(pList{j}) = obj.(pList{j});
            end
        end
    
    end %methods
    
end %class

