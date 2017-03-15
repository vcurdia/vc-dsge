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
        TimeIdx
        TimeStart
        TimeEnd
        T
        SampleStart
        NPreSample
        Model
        Var
        NVar
        Values
        TickLabels
    end
   
    methods
        function obj = Data(fn)
            if nargin>0
                obj.Source = fn;
            end
        end
        
        function Load(obj)
            Raw = importdata(obj.Source);
            Raw.TimeIdx = Raw.textdata(2:end,1)';
            Raw.TimeStart = Raw.TimeIdx{1};
            Raw.TimeEnd = Raw.TimeIdx{end};
            Raw.T = length(Raw.TimeIdx);
            Raw.Var = {Raw.textdata{1,2:end}};
            Raw.NVar = length(Raw.Var);
            Raw.Values = [Raw.data;nan(Raw.T-size(Raw.data,1),Raw.NVar)];
            
            if isempty(obj.TimeStart)
                obj.TimeStart = Raw.TimeStart
            end
            
            if isempty(obj.TimeEnd)
                obj.TimeEnd = Raw.TimeEnd;
            end
            
            obj.TimeIdx = TimeIdxCreate(obj.TimeStart,obj.TimeEnd);
            [tfDates,idxDates] = ismember(obj.TimeIdx,Raw.TimeIdx);
            if ~all(tfDates)
                error('Data require time periods outside data file!')
            end
            obj.T = length(obj.TimeIdx);
            
            if isempty(obj.SampleStart)
                obj.SampleStart = obj.TimeStart;
            end
            obj.NPreSample = find(ismember(obj.TimeStart,obj.SampleStart);
            
            if isempty(obj.Model)
                obj.Var = Raw.Var;
            else
                obj.Var = obj.Model.ObsVar.Names;
            end
            obj.NVar = length(obj.Var);
            [tfVar,idxVar] = ismember(obj.Model.ObsVar.Names,Raw.Var);
            if ~all(tfVar)
                error('Observable names do not match csv headers!')
            end
            
            obj.Values = Raw.Values(idxDates,idxVar);

        end
        
    end %methods
    
end %class

