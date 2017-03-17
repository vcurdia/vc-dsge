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
        TimeStart
        TimeEnd
        SampleStart
        Var
    end

    properties (SetAccess=protected)
        Values
        TimeIdx
        T
        NPreSample
        NVar
        Tick
        TickLabels
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
                obj.setTimeIdx(Raw.textdata(2:end,1)')
                obj.setTick
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
        
        function setTimeIdx(obj,tid)
            if nargin==1
                tid = {obj.TimeStart,obj.TimeEnd};
            end
            if length(tid)==2
                tid = TimeIdxCreate(tid{:});
            end
            if ~isempty(obj.TimeIdx)
                obj.Values = obj.Values(ismember(obj.TimeIdx,tid),:);
            end
            obj.TimeIdx = tid;
            obj.T = length(tid);
            if isempty(obj.TimeStart)
                obj.TimeStart = tid{1};
            end
            if isempty(obj.TimeEnd)
                obj.TimeEnd = tid{end};
            end
            if isempty(obj.SampleStart) ...
                    || ~ismember(obj.SampleStart,obj.TimeIdx)
                obj.SampleStart = obj.TimeStart;
            end
        end
        
        function set.TimeStart(obj,t)
            if ~ismember(t,obj.TimeIdx)
                error('Requested TimeStart out of data scope.')
            end
            obj.TimeStart = t;
            if ~isempty(obj.TimeEnd)
                obj.setTimeIdx
            end
        end
        
        function set.TimeEnd(obj,t)
            if ~ismember(t,obj.TimeIdx)
                error('Requested TimeEnd out of data scope.')
            end
            obj.TimeEnd = t;
            if ~isempty(obj.TimeStart)
                obj.setTimeIdx 
            end
        end
        
        function set.SampleStart(obj,t)
            if ~ismember(t,obj.TimeIdx)
                error('Requested SampleStart out of data scope.')
            end
            obj.SampleStart = t;
            obj.NPreSample = find(ismember(obj.TimeIdx,obj.SampleStart))-1;
        end
        
        function setTickLabels(obj,tDates)
            tDates = obj.TimeIdx(ismember(obj.TimeIdx,tDates));
            obj.TickLabels = tDates;
            if isempty(obj.Tick) 
                obj.setTick
            end
            if ~all(ismember(tDates,obj.TimeIdx(obj.Tick)))
                obj.setTick(obj.TickLabels)
            end
            tDates = obj.TimeIdx(obj.Tick);
            tDates(~ismember(tDates,obj.TickLabels)) = {''};
            obj.TickLabels = tDates;
        end
        
        function setTick(obj,idx)
            if nargin<2
                idx = obj.findfirstq4:4:obj.T;
            end
            if iscell(idx)
                idx = find(ismember(obj.TimeIdx,idx));
            end
            idx = idx(idx>0);
            idx = idx(idx<obj.T);
            obj.Tick = idx;
            if isempty(obj.TickLabels)
                tStep = ceil(obj.T/4);
                while mod(tStep,4), tStep = tStep-1; end
                obj.setTickLabels(obj.TimeIdx(obj.findfirstq4:tStep:obj.T))
            end
            TickLabels = obj.TimeIdx(idx);
            TickLabels(~ismember(TickLabels,obj.TickLabels)) = {''};
            obj.TickLabels = TickLabels;
        end
        
        function t = findfirstq4(obj)
            for t=1:obj.T
                if strcmp(obj.TimeIdx{t}(end),int2str(4))
                    break
                end
            end
        end
        
        function new = copy(obj)
            new = DSGE.Data;
            % Copy all non-hidden properties.
            pList = properties(obj);
            for i = 1:length(pList)
                new.(pList{i}) = obj.(pList{i});
            end
        end
    
    end %methods
    
end %class


