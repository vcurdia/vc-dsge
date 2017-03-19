function out = SumStats(data,varargin)

% SumStats
%
% Descriptive statistics
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%% Default Options
op.Percentiles = [0.01, 0.025, 0.05, 0.15, 0.25, ...
                       0.75, 0.85, 0.95, 0.975, 0.99];

%% Update options
op = UpdateOptions(op,varargin{:});

%% create stats
out = struct;
out.Mean = mean(data,2);
out.SD = std(data,0,2);
out.Min = min(data,2);
out.Max = max(data,2);
out.Median = prctile(data,50,2);
for jPrc=1:length(op.Percentiles)
    out.(op.Percentiles{jPrc}) = ...
        prctile(data,100*op.Percentiles(jPrc),2);
end

