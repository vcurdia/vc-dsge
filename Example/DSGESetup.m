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

% diary('MyDSGE.log')
% diary on

% parpool(2)

%% -------------------------------------------------------------------

%% Initiate DSGE

OpenDSGE('MyDSGE')


%% Setup the model

Model.ObsVar = {'DGDP';'PI';'FFR'};

Model.StateVar = {...
    % Regular variables
    'xtil';'YA';'pitil';'pi';'ir';'r';
    'xe';'re';'YAe';
    'delta';'gamma';'u';
    % a couple of artificial variables
    'YAL';'YAeL';
    };

Model.ShockVar = {'edelta';'egamma';'eu';'ei'}; 

Model.AuxParam = {...
    'beta','0.99';
    'gamma','gammaa/400';
    'r','ra/400';
    'phigammatil','exp(gamma)/(exp(gamma)-beta*eta)';
    'etagammatil','exp(gamma)/(exp(gamma)-eta)';
    'phigamma','phigammatil*etagammatil';
    'etagamma','eta/exp(gamma)';
    };

Model.ObsEq = {...
    'gammaa*one+400*(YA_t-YAL_t+gamma_t) - DGDP_t';
    'pistar*one+400*pi_t - PI_t';
    '(ra+pistar)*one+400*ir_t - FFR_t';
    };

Model.StateEq = {...
    % IS Block
    'xtil_tF-phigamma^(-1)*(ir_t-pi_tF-re_t)-xtil_t';
    ['(xe_t-etagamma*(YAL_t-YAeL_t))-beta*etagamma*(xe_tF-etagamma*xe_t)',...
     '-xtil_t'];
    'ir_t-pi_tF - r_t';
    % Efficient Rates
    'YA_t-YAe_t-xe_t';
    'YAe_tF-omega^(-1)*(gamma_tF-re_t+delta_tF)-YAe_t';
    ['-phigamma*(YAe_t-etagamma*(YAeL_t-gamma_t)',...
     '-beta*etagamma*(YAe_tF+gamma_tF-etagamma*YAe_t))',...
     '+beta*etagamma/(1-beta*etagamma)*delta_tF-omega*YAe_t'];
    % PC Block
    'beta*pitil_tF+xi*(omega*xe_t+phigamma*xtil_t)+u_t-pitil_t';
    'pi_t-zeta*pi_tL-pitil_t';
    % Policy Rule
    'rho*ir_tL+(1-rho)*(phipi*pi_t+phix/4*xe_t)+sigmai/400*ei_t-ir_t';
    % Shocks
    'rhodelta*delta_tL+sigmadelta/400*edelta_t-delta_t';
    'rhogamma*gamma_tL+sigmagamma/400*egamma_t-gamma_t';
    'rhou*u_tL+sigmau/400*eu_t-u_t';
    % Auxiliary equations
    'YAL_t-YA_tL'; 
    'YAeL_t-YAe_tL';
    };





%% -------------------------------------------------------------------

%% Finish up

CloseDSGE
fprintf('\n%s\n\n',vctoc)

% delete(gcp)
% diary off
% exit

%% -------------------------------------------------------------------
