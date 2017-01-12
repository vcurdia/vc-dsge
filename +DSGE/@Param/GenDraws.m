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
        nDraws = 1;
    end
    
    if strcmp(Dist,'PriorDraws')
        xd = obj.Prior.Draw(nDraws);
    elseif strcmp(Dist,'PostDraws')
        xd = obj.Post.Draw(nDraws);
%         load(obj.FileName.MCMCDrawsRedux,'xd')
%         xd = xd(:,randi(size(xd,2),1,nDraws));
    else
        if ~isfield(obj,Dist)
            error('Did not recognize distribution to use.\n');
            return
        end
        xd = repmat(obj.(Dist),1,nDraws);
    end

end

%% -------------------------------------------------------------------
