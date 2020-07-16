function run(obj,varargin)

% run
% 
% Run MCMC chains
%
% see also:
% DSGE.MCMC
%
% Created: July 16, 2020
% Copyright (C) 2020 Vasco Curdia

    %% allow updated options to mcmc to be passed as arg to run
    updateoptions(obj,varargin{:});
    
    %% check chains
    obj.checkchains
    
    %% initiate chains
    if size(obj.x0,2)<obj.nchains
        obj.initchains
    end
    

end
