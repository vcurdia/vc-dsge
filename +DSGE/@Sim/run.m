function run(obj,dist,ndraws)

% run
% 
% sun model simulations
% 
% See also:
% DSGE.Sim
%
% ............................................................................
% 
% Created: May 2, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Preamble
if nargin>1, obj.Dist = dist; end
if nargin>2, obj.NDraws = ndraws; end

%% Prepare parameter draws
if strcmp(obj.Dist,'PriorDraws')
    xd = obj.Prior.draw(obj.NDraws);
elseif ismember(obj.Dist,{'PriorMean','PriorMode','PriorMedian'})
    xd = repmat(obj.Prior.(obj.Dist{6:end}),1,obj.NDraws);
elseif strcmp(obj.Dist,'PostDraws')
    xd = obj.Post.draw(obj.NDraws);
elseif ismember(obj.Dist,{'PostMean','PostMode','PostMedian'})
    xd = repmat(obj.Post.(obj.Dist{5:end}),1,obj.NDraws);
else
    xd = repmat(obj.Model.Param.Values,1,obj.NDraws);
end
simName = [obj.Name,'_',obj.Dist];

%% run sims
fprintf('\n*** Running simulations\n')
fprintf('%s\n',simName)
obj.TimeTracker.start(simName)


obj.irf(xd,op)
obj.vd(xd,op)
obj.states(data,xd,op)
obj.sd(data,xd,op)

%% Finish up
close all
obj.TimeTracker.stop(simName)

