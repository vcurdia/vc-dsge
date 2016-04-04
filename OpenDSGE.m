function dsge = OpenDSGE(Spec,SpecPath,SpecPathBase)

% OpenDSGE
%
% Usage:
%   dsge = CloseDSGE(Spec)
%   dsge = CloseDSGE(Spec, SpecPath)
%   dsge = CloseDSGE(Spec, SpecPath, SpecPathBase)
%
% Changes current folder to that of SpecPath and Initializes or loads DSGE as a
% structure.
%
% Note: the structure with the DSGE in the workspace does not have to be called
%       "dsge".
%
% See also:
% RunMyDSGE, CloseDSGE
%
% ...........................................................................
%
% Created: January 22, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

if ~exist('SpecPath','var')
    SpecPath = [Spec,'/'];
    SpecPathBase = '../';
end
if ~exist('SpecPathBase','var')
    if strcmp(SpecPath,'./')
        SpecPathBase = './';
    else
        SpecPathBase = pwd;
    end
end

if ~isempty(SpecPath) && ~isdir(SpecPath)
    mkdir(SpecPath)
end
cd(SpecPath)

if exist([Spec,'.mat'],'file')
    fprintf('Opening spec %s\n',Spec);
    dsge = load(Spec);
else
    fprintf('Initiating DSGE: %s\n',Spec);
    dsge.Spec = Spec;
    dsge.SpecPath = SpecPath;
    dsge.SpecPathBase = SpecPathBase;
    dsge.FileName = struct;
    dsge.Report = struct;
    dsge.Options = struct;
    dsge.Status = struct;
    dsge.TimeElapsed = struct;
    save(Spec,'-struct','dsge')
end

