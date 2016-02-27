function s = OpenDSGE(Spec,SpecPath,SpecPathBase)

% OpenDSGE
%
% Usage:
%   s = CloseDSGE(Spec)
%   s = CloseDSGE(Spec, SpecPath)
%   s = CloseDSGE(Spec, SpecPath, SpecPathBase)
%
% Changes current folder to that of SpecPath and Initializes or loads DSGE as a
% structure.
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
    s = load(Spec);
else
    fprintf('Initiating DSGE: %s\n',Spec);
    s.Spec = Spec;
    s.SpecPath = SpecPath;
    s.SpecPathBase = SpecPathBase;
    s.FileName = struct;
    s.Report = struct;
    s.Options = struct;
    s.Status = struct;
    s.TimeElapsed = struct;
    save(Spec,'-struct','s')
end

