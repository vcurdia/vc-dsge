function sim(obj,varargin)

% sim
% 
% Simulate DSGE model
% 
% See also:
% DSGE.Model, setupdsge
%
% .............................................................................
% 
% Created: May 1, 2017 by Vasco Curdia
% Copyright 2017-2018 by Vasco Curdia


%% Default Options
op.FNSuffix = '';
% op.List = {'IRF','VD','States','SD'};
op.List = {'IRF','States'};
op.Dist = '';
op.NDraws = [];
op.Draws = [];
op = updateoptions(op,varargin{:});


%% prepare draws to use
if ~isempty(op.Draws)
    xd = op.Draws;
else
    if isempty(op.NDraws) 
        if ismember(op.Dist,{'PriorDraws','PostDraws'})
            op.NDraws = 1000;
        else
            op.NDraws = 1;
        end
    end
    isValues = 0;
    if strcmp(op.Dist,'PriorDraws')
        xd = obj.priordraw(op.NDraws);
    elseif ismember(op.Dist,{'PriorMean','PriorMode','PriorMedian'})
        xd = repmat(obj.Prior.(op.Dist(6:end)),1,op.NDraws);
    elseif strcmp(op.Dist,'PostDraws')
        xd = obj.postdraw(op.NDraws);
    elseif ismember(op.Dist,{'PostMean','PostMode','PostMedian'})
        xd = repmat(obj.Post.(op.Dist(5:end)),1,op.NDraws);
    else
        isValues = 1;
        xd = repmat(obj.Param.Values,1,op.NDraws);
    end
    if isempty(op.FNSuffix) && ~isValues
        op.FNSuffix = ['-',lower(op.Dist)];
    end
end


%% run sim
if ismember('IRF',op.List), obj.irf(xd,op); end
if ismember('VD',op.List), obj.vd(xd,op), end
if ismember('States',op.List), obj.states(xd,op), end
if ismember('SD',op.List), obj.sd(xd,op), end


end

