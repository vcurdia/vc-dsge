function SaveDSGE(dsge)

% SaveDSGE
%
% Usage:
%   SaveDSGE(dsge)
%
% Saves DSGE (but does not leave dsge folder).
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% RunMyDSGE, OpenDSGE, CloseDSGE
%
% ...........................................................................
%
% Created: October 20, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%fprintf('\nSaving workspace and exiting spec folder\n')
save(dsge.Spec,'-struct','dsge')


