function TableBreaks = SetTableBreaks(obj,n)

% SetTableBreaks
%
% Set index of table breaks
%
% ...........................................................................
% 
% Created: January 16, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

TableBreaks = obj.TableMaxRows:obj.TableMaxRows:n;
if ~ismember(n,TableBreaks)
    TableBreaks(end+1) = n; 
end
