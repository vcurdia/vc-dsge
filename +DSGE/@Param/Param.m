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
        Prior
        Post
        IsCalibrated
        EstimateIdx
        Percentiles = [0.01, 0.025, 0.05, 0.15, 0.25, ...
                       0.75, 0.85, 0.95, 0.975, 0.99];
        PriorNDraws = 1000;
        TablePrecision = 3;
        TableMaxRows = 35;
        TableMoveLeft = 1; 
        TableLines = [];
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
                    obj.Prior.Dist = p(:,2);
                    obj.Prior.Mean = [p{:,3}]';
                    for j=1:obj.N
                        if isempty(p{j,4})
                            p{j,4} = 0;
                        end
                    end
                    obj.Prior.SD = [p{:,4}]';
                    obj.Values = obj.Prior.Mean;
                    obj.EstimateIdx = ~ismember(obj.Prior.Dist,{'C'});
                    obj.IsCalibrated = ~any(obj.EstimateIdx);
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

