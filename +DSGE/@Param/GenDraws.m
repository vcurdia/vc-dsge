function xd = GenDraws(obj,Dist,nDraws)

% GenDraws
% 
% Generates a sample of draws for the parameter vector.
% 
% ...........................................................................
% 
% Created: January 12, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia

%----------------------------------------------------------------------------

xd = [];
if nargin<1 || isempty(Dist)
    Dist = 'Values';
end
if nargin<2 || isempty(nDraws)
    if ismember(Dist,{'PriorDraws','PostDraws'})
        nDraws = 1000; 
    else
        nDraws = 1; 
    end
end

if strcmp(Dist,'PriorDraws')
    xd = obj.DrawPrior(nDraws);
elseif strcmp(Dist,'PostDraws')
    xd = obj.Post.Draw(nDraws);
%         load(obj.FileName.MCMCDrawsRedux,'xd')
%         xd = xd(:,randi(size(xd,2),1,nDraws));
else
    if ~isprop(obj,Dist)
        error('Did not recognize distribution to use.');
        return
    end
    if isempty(obj.(Dist))
        error('%s is empty.',Dist);
        return
    end
    xd = repmat(obj.(Dist),1,nDraws);
end

