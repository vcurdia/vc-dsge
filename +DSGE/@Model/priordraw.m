function xd = priordraw(obj,nDraws)

% drawprior
%
% Pick parameters from the prior distribution sample 
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
if isempty(obj.Prior.Sample.FileName)
    obj.makepriorsample
end
draws = load(obj.Prior.Sample.FileName);
xd = draws.Param(:,randi(draws.N,1,nDraws));
