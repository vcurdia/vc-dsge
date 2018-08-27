function p = priorlpdf(obj,x)

% priorlpdf
%
% evaluate the prior log-pdf for a set of parameters 
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% Copyright 2017-2018 by Vasco Curdia


p = log(obj.priorpdf(x));
