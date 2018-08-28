% setupdsge
%
% Example of how to setup a DSGE using the VC-DSGE toolbox.
%
%
% See also:
% estimatedsge, DSGE.Model, DSGE.Data
%
% Created: January 21, 2016 by Vasco Curdia
% Copyright 2016-2018 Vasco Curdia


%% Preamble
clear all
set(groot,'defaultTextInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
setpath

%% Initiate time tracker
tt = TimeTracker;
tt.start('Setup')

%% Setup the model
model = DSGE.Model('MyDSGE');
mkdir(model.Name)
cd(model.Name)

% % example for calibrated model
% model.Param = {...
%     'beta', 0.99,'$\beta$';
%     'omega', 1,'$\omega$';
%     'xi', 0.1, '$\xi$';
%     'eta', 0.6, '$\eta$';
%     'zeta', 0.6,'$\zeta$';
%     'rho', 0.7, '$\rho$';
%     'phipi', 1.5, '$\phi_\pi$';
%     'phix', 0.5, '$\phi_x$';
%     'pistar', 2, '$\pi^*$';
%     'ra', 2,'$r^a$';
%     'gammaa', 3,'$400\gamma$';
%     'rhodelta', 0.5, '$\rho_\delta$';
%     'rhogamma', 0.5, '$\rho_\gamma$';
%     'rhou', 0.5, '$\rho_u$';
%     'sigmadelta', 0.5,'$\sigma_\delta$';
%     'sigmagamma', 0.5,'$\sigma_\gamma$';
%     'sigmau', 0.5, '$\sigma_u$';
%     'sigmai', 0.5, '$\sigma_i$';
%                  };

% example for model w/ prior, to be estimated
model.Param = {...
    'beta', 'C', 0.99, [], '$\beta$';
    'omega', 'G', 1, 0.2, '$\omega$';
    'xi', 'G', 0.1, 0.05, '$\xi$';
    'eta', 'B', 0.6, 0.2, '$\eta$';
    'zeta', 'B', 0.6, 0.2, '$\zeta$';
    'rho', 'B', 0.7, 0.15, '$\rho$';
    'phipi', 'N', 1.5, 0.25, '$\phi_\pi$';
    'phix', 'N', 0.5, 0.2, '$\phi_x$';
    'pistar', 'N', 2, 1, '$\pi^*$';
    'ra', 'N', 2, 1, '$r^a$';
    'gammaa', 'N', 3, 0.35, '$400\gamma$';
    'rhodelta', 'B', 0.5, 0.2, '$\rho_\delta$';
    'rhogamma', 'B', 0.5, 0.2, '$\rho_\gamma$';
    'rhou', 'B', 0.5, 0.2, '$\rho_u$';
    'sigmadelta', 'IG1', 0.5, 2, '$\sigma_\delta$';
    'sigmagamma', 'IG1', 0.5, 2, '$\sigma_\gamma$';
    'sigmau', 'IG1', 0.5, 2, '$\sigma_u$';
    'sigmai', 'IG1', 0.5, 2, '$\sigma_i$';
                 };

% Uncomment the following lines to show how Param.NumSolve works: 
model.NumSolveParam = {...
    'rA',1,'$r^A$';
    'rB',1,'$r^B$';
                   };
model.NumSolveEq = {...
    '(rA+rB)/2-r';
    'rA+0.5/400-rB';
               };

model.CompositeParam = {...
    'gamma','gammaa/400','$\gamma$';
    'r','ra/400','$r$';
    'phigammatil','exp(gamma)/(exp(gamma)-beta*eta)','$\tilde{\phi}_\gamma$';
    'etagammatil','exp(gamma)/(exp(gamma)-eta)','$\tilde{\eta}_\gamma$';
    'phigamma','phigammatil*etagammatil','$\phi_\gamma$';
    'etagamma','eta/exp(gamma)','$\eta_\gamma$';
                   };

model.ObsVar = {...
    'DGDP', 'GDP Growth';
    'PI', 'Inflation';
    'FFR', 'FFR';
                };

model.StateVar = {...
    % Regular variables
    'xtil', '$\tilde{x}$';
    'YA', '$Y_A$';
    'pitil', '$\tilde{\pi}$';
    'pi', '$\pi$';
    'ir', '$i$';
    'xe', '$x^e$';
    're', '$r^e$';
    'YAe', '$Y_A^e$';
    'delta', '$\delta$';
    'gamma', '$\gamma$';
    'u', '$u$';
    % a couple of artificial variables
    'YAL', '$Y_{A,t-1}$';
    'YAeL', '$Y_{A,t-1}^e$';
                  };

model.ShockVar = {...
    'edelta', '$\varepsilon_\delta$';
    'egamma', '$\varepsilon_\gamma$';
    'eu', '$\varepsilon_u$';
    'ei', '$\varepsilon_i$';
                  }; 

model.AuxVar = {'r','ir_t-pi_tF','$r$'};

model.ObsEq = {...
    'gammaa+400*(YA_t-YAL_t+gamma_t) - DGDP_t';
    'pistar+400*pi_t - PI_t';
    'ra+pistar+400*ir_t - FFR_t';
               };

model.StateEq = {...
    % IS Block
    'xtil_tF-phigamma^(-1)*(ir_t-pi_tF-re_t)-xtil_t';
    ['(xe_t-etagamma*(YAL_t-YAeL_t))-beta*etagamma*(xe_tF-etagamma*xe_t)',...
     '-xtil_t'];
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

%% Generate Mats
% model.genmats([basePath,'Mats/'])
model.genmats
mats = model.mats(model.Param.Values);

%% Analyze Prior parameters
model.initializeprior

%% Data
% data = DSGE.Data([basePath,'Data/Data_1987q3_2009q3.csv']);
model.Data = DSGE.Data('../Data_1987q3_2009q3.csv');
model.Data.Var = model.ObsVar.Names;
op.Sim.Data = model.Data;
op.Sim.Tick.Labels = {'1990q1','1995q1','2000q1','2005q1'};

%% Create posterior
model.initializepost

%% Simulate model at prior
% model.sim(prior.draw(1000),op.Sim,'FNSuffix','_PriorDraws')

%% Finish up
tt.stop('Setup')
save(model.Name)
cd ..
