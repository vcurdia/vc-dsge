% DSGESetup
%
% This file gives an example of how to use the VC_DSGE package
%
% Convention: x_t refers to x(t)
%             x_tF refers to x(t+1)
%             x_tL refers to x(t-1)
%             x_ss refers to steady state of x(t)
%
%   
% See also:
%
% ...........................................................................
%
% Created: January 21, 2016 by Vasco Curdia
% 
% Copyright (C) 2016 Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble
clear all
tic

%% Settings
Spec = 'MyDSGE';


%% Initiate parpool
% parpool(2)

%% -------------------------------------------------------------------

%% create folder if needed
if ~isdir(Spec)
    mkdir(Spec)
end
cd(Spec)
fprintf('\nSpec: %s\n',Spec);


%% -------------------------------------------------------------------

%% Close parpool if needed
% delete(gcp)

%% Save workspace
save('Output')
cd ..
    
%% elapsed time
fprintf('\n%s\n\n',vctoc)

%% exit if in server
% exit

%% -------------------------------------------------------------------
