function setpath

% setpath
%
% Set path to codes needed.
%
% Created: August 18, 2016 by Vasco Curdia

%% set path
pathBase = fullfile(getenv('HOME'),'matlab');
pathList = {...
    'vc-dsge',...
    'vc-tools',...
    'sims/gensys',...
    'sims/kf',...
    'sims/optimize',...
    };
pathAdd = fullfile(pathBase,pathList);
addpath(pathAdd{:})

%% set text interpreter to latex
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
