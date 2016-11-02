function obj = CreateDSGE(Spec,SpecPath)

% createdsge
%
% Usage:
%   obj = CreateDSGE(Spec)
%   obj = CreateDSGE(Spec, SpecPath)
%
% Creates subfolder under SpecPath with name of Spec and initializes the DSGE 
% as a structure.
%
% Inputs:
%   Spec: string with spec name, used for folder name to hold all the spec data 
%         and spec-specific functions
%   SpecPath: base path where spec folder will be located. Needs to end in "/".
% 
% Note: the structure with the DSGE in the workspace does not have to be called
%       "obj".
%
% See also:
% SetupMyDSGE, OpenDSGE, SaveDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

if exist('SpecPath','var')
    if ~isdir(SpecPath), mkdir(SpecPath), end
    cd(SpecPath)
end
mkdir(Spec)
cd(Spec)

fprintf('Initiating DSGE: %s\n',Spec);
obj.Spec = Spec;
obj.PlotDir = struct;
obj.Options = struct;
obj.TimeElapsed = struct;
obj.Param = struct;
obj.Var = struct;
obj.Eq = struct;

SaveDSGE(obj)


