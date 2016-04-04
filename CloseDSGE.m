function dsge = CloseDSGE(dsge)

% CloseDSGE
%
% Usage:
%   dsge = CloseDSGE(dsge)
%   dsge = CloseDSGE(dsge,BaseFolder)
%
% Saves DSGE and exists Spec folder.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% RunMyDSGE, OpenDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%fprintf('\nSaving workspace and exiting spec folder\n')
save(dsge.Spec,'-struct','dsge')
cd(dsge.SpecPathBase)

