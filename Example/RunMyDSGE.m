% RunMyDSGE
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
% OpenDSGE, CloseDSGE, PrepModel
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

BaseFolder = pwd;

% diary('MyDSGE.log')
% diary on

% parpool(2)

%% -------------------------------------------------------------------

%% Initiate DSGE
s = OpenDSGE('MyDSGE');


%% Setup the model

Model.Param = {...
    'omega', 'G', 1, 0.2,'\omega';
    'xi', 'G', 0.1, 0.05,'\xi';
    'eta', 'B', 0.6, 0.2,'\eta';
    'zeta', 'B', 0.6, 0.2,'\zeta';
    'rho', 'B', 0.7, 0.15,'\rho';
    'phipi', 'N', 1.5, 0.25,'\phi_\pi';
    'phix', 'N', 0.5, 0.2,'\phi_x';
    'pistar', 'N', 2, 1,'\pi^*';
    'ra', 'N', 2, 1,'r^a';
    'gammaa', 'N', 3, .35,'\gamma^a';
    'rhodelta', 'B', 0.5, 0.2,'\rho_\delta';
    'rhogamma', 'B', 0.5, 0.2,'\rho_\gamma';
    'rhou', 'B', 0.5, 0.2,'\rho_u';
    'sigmadelta', 'IG1', 0.5, 2,'\sigma_\delta';
    'sigmagamma', 'IG1', 0.5, 2,'\sigma_\gamma';
    'sigmau', 'IG1', 0.5, 2,'\sigma_u';
    'sigmai', 'IG1', 0.5, 2,'\sigma_i';
    };

Model.AuxParam = {...
    'beta','0.99','\beta';
    'gamma','gammaa/400','\gamma';
    'r','ra/400','r';
    'phigammatil','exp(gamma)/(exp(gamma)-beta*eta)','phigammatil';
    'etagammatil','exp(gamma)/(exp(gamma)-eta)','etagammatil';
    'phigamma','phigammatil*etagammatil','\phi_\gamma';
    'etagamma','eta/exp(gamma)','\eta_\gamma';
    };

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

s = PrepModel(s,Model);



%% -------------------------------------------------------------------

%% Finish up
s = CloseDSGE(s,BaseFolder);
fprintf('\n%s\n\n',vctoc)

% delete(gcp)
% diary off
% exit

%% -------------------------------------------------------------------
