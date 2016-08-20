function dsge = CreateDSGE(Spec,SpecPath,SpecPathBase)

% CreateDSGE
%
% Usage:
%   dsge = CloseDSGE(Spec)
%   dsge = CloseDSGE(Spec, SpecPath)
%   dsge = CloseDSGE(Spec, SpecPath, SpecPathBase)
%
% Changes current folder to that of SpecPath and Initializes DSGE as a
% structure.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% RunMyDSGE, OpenDSGE, CloseDSGE
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
if ~exist('SpecPathBase','var')
    if strcmp(SpecPath,'./')
        SpecPathBase = './';
    else
        SpecPathBase = '../';
    end
end

if ~isempty(SpecPath) && ~isdir(SpecPath)
    mkdir(SpecPath)
end
cd(SpecPath)

fprintf('Initiating DSGE: %s\n',Spec);
dsge.Spec = Spec;
dsge.SpecPath = SpecPath;
dsge.SpecPathBase = SpecPathBase;
dsge.FileName = struct;
dsge.Report = struct;
dsge.PlotDir = struct;
dsge.Options = struct;
dsge.Status = struct;
dsge.TimeElapsed = struct;
dsge.Param = cell(0,4);
dsge.AuxParam = cell(0,2);
dsge.ObsVar = cell(0,1);
dsge.StateVar = cell(0,1);
dsge.ShockVar = cell(0,1);
dsge.AuxVar = cell(0,2);
dsge.ObsEq = cell(0,1);
dsge.StateEq = cell(0,1);
dsge.AuxEq = cell(0,1);

save(Spec,'-struct','dsge')


