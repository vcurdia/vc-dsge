function initializeposterior(obj)

% initializeposterior
%
% initialize posterior distribution structure
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: August 27, 2018 by Vasco Curdia
% Copyright 2018 by Vasco Curdia


fprintf('Initializing posterior\n')

obj.Posterior.EstimateIdx = ~ismember(obj.Prior.Dist,{'C'});
obj.Posterior.NEstimate = sum(obj.Posterior.EstimateIdx);
obj.Posterior.Mode = obj.Param.Values;
obj.Posterior.ModeLPDF = obj.posteriorlpdf(obj.Posterior.Mode);
fprintf('Posterior log-pdf using Param.Values is %0.4f.\n',...
        obj.Posterior.ModeLPDF);
obj.Posterior.Mean = obj.Prior.Mean;
obj.Posterior.SD = obj.Prior.SD;
obj.Posterior.Var = diag(obj.Prior.SD.^2);
obj.Posterior.Median = obj.Prior.Median;
obj.Posterior.Prc05 = obj.Prior.Prc05;
obj.Posterior.Prc95 = obj.Prior.Prc95;
