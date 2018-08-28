function xd = priordraw(obj,nDraws)

% drawprior
%
% draw random parameter vector(s) from the prior distribution 
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% Copyright 2017-2018 by Vasco Curdia


if nargin<2 || isempty(nDraws)
    nDraws = 1; 
end
xd = nan(obj.Param.N,nDraws);
for j=1:obj.Param.N
    xd(j,:) = obj.Prior.RndCmd{j}(nDraws);
end
