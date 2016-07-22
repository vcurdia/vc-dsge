function dsge = OpenDSGE(Spec)

% OpenDSGE
%
% Usage:
%   dsge = CloseDSGE(Spec)
%   dsge = CloseDSGE(Spec, SpecPath)
%
% Changes current folder to that of SpecPath and loads DSGE as a structure.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% RunMyDSGE, CreateDSGE, CloseDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

if ~exist('SpecPath','var')
    SpecPath = [Spec,'/'];
end
cd(SpecPath)
fprintf('Opening %s\n',Spec);
dsge = load(Spec);

