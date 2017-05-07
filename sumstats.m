function dStats = sumstats(data,Percentiles)

% sumstats
%
% Generate descriptive statistics
%
% Usage:
%   dStats = sumstats(data)
%   dStats = sumstats(data,Percentiles)
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%% Default Options
if nargin<2
    Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
end

%% create stats
dStats = struct;
dStats.Mean = mean(data,2);
dStats.SD = std(data,0,2);
dStats.Min = min(data,[],2);
dStats.Max = max(data,[],2);
dStats.Median = prctile(data,50,2);
for jPrc=1:length(Percentiles)
    dStats.(sprintf('Prc%02.0f',100*Percentiles(jPrc))) = ...
        prctile(data,100*Percentiles(jPrc),2);
end

