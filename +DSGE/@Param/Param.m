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
        PriorMode
        PriorMedian
        PostMean
        PostSD
        PostMode
        PostMedian
        Prior
        Post
        IsCalibrated
        EstimateIdx
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
                    obj.IsCalibrated = 1;
                    obj.Values = [p{:,2}]';
                else
                    obj.PriorDist = p(:,2);
                    obj.EstimateIdx = ~ismember(obj.PriorDist,{'C'});
                    obj.IsCalibrated = ~any(obj.EstimateIdx);
                    obj.PriorMean = p(:,3);
                    obj.PriorSD = p(:,4);
                    obj.Values = obj.PriorMean;
                end
            end
        end
        
        function Show(obj,PropName)
            p = struct;
            if nargin==0 || isempty(PropName)
                PropName = 'Values'; 
            end
            for j=1:obj.N
                p.(obj.Names{j}) = obj.(PropName)(j);
            end
            fprintf('Parameter %s:\n',PropName)
            disp(p)
        end
        
        function SetValues(obj,ParNames,ParValues)
            if ischar(ParNames), ParNames = {ParNames}; end
            np = length(ParNames);
            nv = length(ParValues);
            if np>0 && np==nv
                [tf,idxp] = ismember(ParNames,obj.Names);
                obj.Values(idxp(tf)) = ParValues(tf);
                if ~all(tf)
                    fprintf('Invalid parameter names ignored:\n')
                    fprintf('  %s\n',ParNames{~tf})
                end
            else
                error(['List of parameter names is empty or has different ' ...
                       'length from list of values.'])
            end
        end
        
        
    end %methods
    
end %class

