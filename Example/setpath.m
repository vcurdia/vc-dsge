function setpath

% setpath
%
% Set path to codes needed.
%
% Created: August 18, 2016 by Vasco Curdia
% Copyright 2016-2018 by Vasco Curdia

%% set path
pathBase = '../../../Matlab/';
pathList = {...
    'VC-DSGE',...
    'VC-Tools',...
    'Sims-Gensys',...
    'Sims-KF',...
    'Sims-Optimize',...
    };
for j=1:length(pathList)
    pathAdd{j} = [pathBase,pathList{j}];
end
addpath(pathAdd{:})

%% set text interpreter to latex
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
