function dStats = genstats(data,varargin)

% genstats
%
% Generate descriptive statistics
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
dStats = struct;
dStats.Mean = mean(data,2);
dStats.SD = std(data,0,2);
dStats.Min = min(data,2);
dStats.Max = max(data,2);
dStats.Median = prctile(data,50,2);
for jPrc=1:length(op.Percentiles)
    dStats.(op.Percentiles{jPrc}) = ...
        prctile(data,100*op.Percentiles(jPrc),2);
end

