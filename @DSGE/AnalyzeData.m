function obj = AnalyzeData(obj)

% Analyze the data used for the estimation
%
% See also:
% DSGE, SetupMyDSGE
%
% ...........................................................................
%
% Created: December 13, 2016 by Vasco Curdia
%
% Copyright 2016 by Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

fprintf('\n*** Analyzing Data\n')

%% Load data
Data = importdata(obj.FileName.Data);
Data.Var = Data.textdata(1,2:end);
Data.TimeIdx = Data.textdata(2:end,1)';
[tfDates,idxDates] = ismember(obj.TimeIdx,Data.TimeIdx);
if ~all(tfDates)
    error('DateLabels require time periods outside data file!')
end
[tfVar,idxVar] = ismember(obj.ObsVar.Names,Data.Var);
if ~all(tfVar)
    error('Observable names do not match csv headers!')
end
obj.Data = Data.data(idxDates,idxVar);

%% -------------------------------------------------------------------

end

