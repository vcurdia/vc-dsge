function linkobj(model,prior,post)

% linkobj
%
% link Prior and posterior to model
% 
% ............................................................................
%
% Created: March 23, 2017 by Vasco Curdia
% Copyright 2017 by Vasco Curdia

if nargin>1
    prior.Model = model;
end
if nargin>2
    post.Model = model;
    post.Prior = prior;
end

