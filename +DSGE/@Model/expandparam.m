function xx = expandparam(obj,x)

% expandparam
%
% expand estimated parameter vector(s) to include calibrated parameteres
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: August 27, 2018 by Vasco Curdia
% Copyright 2018 by Vasco Curdia


xx = repmat(obj.Param.Values,1,size(x,2),size(x,3));
xx(obj.Post.EstimateIdx,:,:) = x;
