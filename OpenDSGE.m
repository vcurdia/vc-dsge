function obj = OpenDSGE(Spec,SpecPath)

% OpenDSGE
%
% Usage:
%   obj = OpenDSGE(Spec)
%   obj = OpenDSGE(Spec, SpecPath)
%
% Changes current folder to that of SpecPath and loads DSGE as a structure.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% SetupMyDSGE, CreateDSGE, CloseDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

if exist('SpecPath','var'), cd(SpecPath), end
cd(Spec)
fprintf('Opening %s\n',Spec);
obj = load(Spec);

