classdef Data < handle
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
        Tick
        TickLabels
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
                obj.setTickLabels;
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
                tid = TimeIdxCreate(tid{:});
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
        
        function set.TickLabels(obj,tDates)
            tDates = obj.TimeIdx(ismember(obj.TimeIdx,tDates));
            obj.TickLabels = tDates;
            if isempty(obj.Tick) ...
                    || ~all(ismember(tDates,obj.TimeIdx(obj.Tick)))
                obj.Tick = obj.TickLabels;
            end
            tDates = obj.TimeIdx(obj.Tick);
            tDates(~ismember(tDates,obj.TickLabels)) = {''};
            obj.TickLabels = tDates;
        end
        
        function set.Tick(obj,idx)
            if iscell(idx)
                idx = find(ismember(obj.TimeIdx,idx));
            end
            idx = idx(idx>0);
            idx = idx(idx<obj.T);
            obj.Tick = idx;
            TickLabels = obj.TimeIdx(idx);
            if isempty(obj.TickLabels)...
                    || ~any(ismember(obj.TickLabels,TickLabels))
                obj.TickLabels = TickLabels;
            end
            TickLabels(~ismember(TickLabels,obj.TickLabels)) = {''};
            obj.TickLabels = TickLabels;
        end
        
        function setTickAnnual(obj,q)
            if nargin<2, q = 4; end
            obj.Tick = obj.findQuarter(q):4:obj.T;
        end
        
        function t = findQuarter(obj,q)
            if nargin<2, q = 4; end
            for t=1:obj.T
                if strcmp(obj.TimeIdx{t}(end),int2str(q))
                    break
                end
            end
        end
        
        function setTickLabels(obj,n,q)
            if nargin<2, n = 5; end
            if nargin<3, q = 4; end
            obj.setTickAnnual(q);
            tStep = ceil(obj.T/(n-1));
            while mod(tStep,4), tStep = tStep-1; end
            obj.TickLabels = obj.TimeIdx(obj.findQuarter(q):tStep:obj.T);
        end
        
        function new = Copy(obj)
            new = DSGE.Data;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for j = 1:length(pList)
                new.(pList{j}) = obj.(pList{j});
            end
        end
    
    end %methods
    
end %class


