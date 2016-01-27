function s = OpenDSGE(Spec)

% OpenDSGE
%
% Usage:
%   s = CloseDSGE(Spec)
%
% Changes current folder to that of Spec and Initializes or loads DSGE as a
% structure s.
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

if ~isdir(Spec)
    mkdir(Spec)
end
cd(Spec)

if exist([Spec,'.mat'],'file')
    fprintf('Opening spec %s\n',Spec);
    s = load(Spec);
else
    fprintf('Initiating DSGE: %s\n',Spec);
    s.Spec = Spec;
    s.Status = struct;
    s.Options = struct;
    save(Spec,'-struct','s')
end

