classdef Data < matlab.mixin.Copyable
% DSGE.Data class
% 
% See also:
% SetupMyDSGE
%
% ...........................................................................
%
% Created: March 14, 2017 by Vasco Curdia
% 
% Copyright (C) 2017 Vasco Curdia
    
    properties
        Source
        Values
        TimeIdx
        Var
        SampleStart
    end

    properties (SetAccess=protected)
        TimeStart
        TimeEnd
        T
        NPreSample
        NVar
    end
    
    methods
        function obj = Data(fn)
            if nargin>0
                fprintf('\n*** Loading data from:\n%s\n',fn)
                obj.Source = fn;
                Raw = importdata(obj.Source);
                obj.Var = {Raw.textdata{1,2:end}};
                obj.Values = [...
                    Raw.data;
                    NaN(size(Raw.textdata,1)-1-size(Raw.data,1),obj.NVar)];
                obj.TimeIdx = Raw.textdata(2:end,1)';
            end
        end
        
        function set.Var(obj,Var)
            if ~isempty(obj.Var)
                [tf,idx] = ismember(Var,obj.Var);
                if ~all(tf)
                    error('Variables not found in data.')
                end
                obj.Values = obj.Values(:,idx);
            end
            obj.Var = Var;
            obj.NVar = length(Var);
        end
        
        function set.TimeIdx(obj,tid)
            if length(tid)==2
                tid = timeidx(tid{:});
            end
            if ~isempty(obj.TimeIdx)
                obj.Values = obj.Values(ismember(obj.TimeIdx,tid),:);
            end
            obj.TimeIdx = tid;
            obj.T = length(tid);
            obj.TimeStart = tid{1};
            obj.TimeEnd = tid{end};
            if isempty(obj.SampleStart) ...
                    || ~ismember(obj.SampleStart,obj.TimeIdx)
                obj.SampleStart = obj.TimeStart;
            end
        end
        
        function set.SampleStart(obj,t)
            if ~ismember(t,obj.TimeIdx)
                error('Requested SampleStart out of data scope.')
            end
            obj.SampleStart = t;
            obj.NPreSample = find(ismember(obj.TimeIdx,obj.SampleStart))-1;
        end
        
    end %methods
    
end %class


