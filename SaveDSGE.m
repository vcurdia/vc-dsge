function SaveDSGE(obj)

% SaveDSGE
%
% Usage:
%   SaveDSGE(obj)
%
% Saves DSGE
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "obj".
%
% See also:
% SetupMyDSGE, OpenDSGE, CreateDSGE
%
% ...........................................................................
%
% Created: October 20, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%fprintf('\nSaving workspace and exiting spec folder\n')
%save(obj.Spec,'-struct','obj')
save('DSGEData','-struct','obj')


