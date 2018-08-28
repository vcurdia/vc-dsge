function xd = postdraw(obj,nDraws)

% postdraw
%
% draw from the posterior sample
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: August 27, 2018 by Vasco Curdia
% Copyright 2018 by Vasco Curdia


if nargin<2 || isempty(nDraws)
    nDraws = 1; 
end
if isempty(obj.Post.MCMCSample.FileNameRedux)
    obj.mcmcredux
end
draws = load(obj.Post.MCMCSample.FileNameRedux);
xd = obj.expandparam(draws.Param(:,randi(draws.N,1,nDraws)));
