function s = CloseDSGE(s,BaseFolder)

% CloseDSGE
%
% Usage:
%   s = CloseDSGE(s)
%   s = CloseDSGE(s,BaseFolder)
%
% Saves DSGE and exists Spec folder.
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

fprintf('Saving workspace and exiting DSGE folder\n')
save(s.Spec,'-struct','s')
if nargin>1
    cd(BaseFolder)
else
    cd ..
end

