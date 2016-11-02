function obj = opendsge(Spec,SpecPath)

% OpenDSGE
%
% Usage:
%   obj = opendsge(Spec)
%   obj = opendsge(Spec, SpecPath)
%
% Changes current folder to that of SpecPath and loads DSGE as a structure.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% setupMyDSGE, createdsge, closedsge
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

