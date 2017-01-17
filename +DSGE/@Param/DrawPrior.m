function xd = DrawPrior(obj,nDraws)

% DrawPrior
% 
% Generates a sample of draws for the parameter vector from the prior.
% 
% ...........................................................................
% 
% Created: January 16, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%----------------------------------------------------------------------------

if nargin<1 || isempty(nDraws)
    nDraws = 1; 
end

xd = nan(obj.N,nDraws);

for j=1:obj.N
    xd(j,:) = obj.Prior.RndCmd{j}(nDraws);
end
