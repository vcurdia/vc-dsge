function obj = closedsge(obj)

% closedsge
%
% Usage:
%   obj = closedsge(obj)
%   obj = CloseDSGE(obj,BaseFolder)
%
% Saves obj and exists Spec folder.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "obj".
%
% See also:
% setupMyDSGE, opendsge
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%fprintf('\nSaving workspace and exiting spec folder\n')
save(obj.Spec,'-struct','obj')
cd(obj.SpecPathBase)

