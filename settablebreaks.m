function tableBreaks = settablebreaks(nRows,maxRows)

% settablebreaks
%
% Set index of table breaks.
% 
% Usage:
%   tableBreaks = settablebreaks(nRows)
%   tableBreaks = settablebreaks(nRows,maxRows)
%
% ...........................................................................
% 
% Created: January 16, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

if nargin<2, maxRows = 35; end
tableBreaks = maxRows:maxRows:nRows;
if ~ismember(nRows,tableBreaks)
    tableBreaks(end+1) = nRows; 
end
