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
        SampleStart
        NPreSample
        Var
        Values
        Tick
        TickLabels
    end

    properties (SetAccess=protected)
        TimeIdx
        TimeStart
        TimeEnd
        T
        NVar
    end
    
    methods
        function obj = Data(fn)
            if nargin>0
                obj.Source = fn;
                obj.Load
            end
        end
        
        function Load(obj)
            fprintf('\n*** Loading data\n')
            TimeElapsed = tic;
            Raw = importdata(obj.Source);
            Raw.TimeIdx = Raw.textdata(2:end,1)';
            Raw.TimeStart = Raw.TimeIdx{1};
            Raw.TimeEnd = Raw.TimeIdx{end};
            Raw.T = length(Raw.TimeIdx);
            Raw.Var = {Raw.textdata{1,2:end}};
            Raw.NVar = length(Raw.Var);
            Raw.Values = [Raw.data;nan(Raw.T-size(Raw.data,1),Raw.NVar)];
            
            if isempty(obj.TimeIdx)
                obj.SetTime(Raw.TimeIdx{1,end});
            end
            [tfDates,idxDates] = ismember(obj.TimeIdx,Raw.TimeIdx);
            if ~all(tfDates)
                error('Requested time periods outside data file!')
            end
            obj.T = length(obj.TimeIdx);
            
            if isempty(obj.SampleStart)
                obj.SetPreSample(0)
            end
            obj.NPreSample = ...
                max(0,find(ismember(obj.TimeIdx,obj.SampleStart))-1);
            
            if isempty(obj.Var)
                obj.Var = Raw.Var;
            end
            obj.NVar = length(obj.Var);
            [tfVar,idxVar] = ismember(obj.Var,Raw.Var);
            if ~all(tfVar)
                error('Variable names do not match csv headers!')
            end
            
            obj.Values = Raw.Values(idxDates,idxVar);
            
            if isempty(obj.Tick)
                obj.SetTick
            end
            
            fprintf('DataLoad: '), vctoc(TimeElapsed)
            
        end
        
        function SetTime(obj,TimeStart,TimeEnd)
            obj.TimeStart = TimeStart;
            obj.TimeEnd = TimeEnd;
            obj.TimeIdx = TimeIdxCreate(TimeStart,TimeEnd);
            obj.T = length(obj.TimeIdx);
        end
        
        function SetPreSample(obj,op)
            if ischar(op)
                obj.SampleStart = op;
                obj.NPreSample = ...
                    max(0,find(ismember(obj.TimeIdx,obj.SampleStart))-1);
            else
                obj.SampleStart = obj.TimeIdx(op+1);
                obj.NPreSample = op;
            end
        end
        
        function SetTickLabels(obj,tDates)
            tDates = obj.TimeIdx(ismember(obj.TimeIdx,tDates));
            obj.TickLabels = tDates;
            if isempty(obj.Tick) 
                obj.SetTick
            end
            if ~all(ismember(tDates,obj.TimeIdx(obj.Tick)))
                obj.SetTick(obj.TickLabels)
            end
        end
        
        function SetTick(obj,idx)
            if nargin<2
                idx = obj.FindFirstQ4:4:obj.T;
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
                obj.SetTickLabels(obj.TimeIdx(obj.FindFirstQ4:tStep:obj.T))
            end
            TickLabels = obj.TimeIdx(idx);
            TickLabels(~ismember(TickLabels,obj.TickLabels)) = {''};
            obj.TickLabels = TickLabels;
        end
        
        function t = FindFirstQ4(obj)
            for t=1:obj.T
                if strcmp(obj.TimeIdx{t}(end),int2str(4))
                    break
                end
            end
        end
        
    end %methods
    
end %class


