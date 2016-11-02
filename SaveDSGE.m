function savedsge(obj)

% savedsge
%
% Usage:
%   savedsge(obj)
%
% Saves DSGE
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "obj".
%
% See also:
% setupMyDSGE, opendsge, createdsge
%
% ...........................................................................
%
% Created: October 20, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%fprintf('\nSaving workspace and exiting spec folder\n')
save(obj.Spec,'-struct','obj')


