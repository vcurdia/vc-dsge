function initializepost(obj)

% initializepost
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

obj.Post.EstimateIdx = ~ismember(obj.Prior.Dist,{'C'});
obj.Post.NEstimate = sum(obj.Post.EstimateIdx);
obj.Post.Mode = obj.Param.Values;
obj.Post.ModeLPDF = obj.postlpdf(obj.Post.Mode);
fprintf('Posterior log-pdf using Param.Values is %0.4f.\n',...
        obj.Post.ModeLPDF);
obj.Post.Mean = obj.Prior.Mean;
obj.Post.SD = obj.Prior.SD;
obj.Post.Var = diag(obj.Prior.SD.^2);
obj.Post.Median = obj.Prior.Median;
obj.Post.Prc05 = obj.Prior.Prc05;
obj.Post.Prc95 = obj.Prior.Prc95;
