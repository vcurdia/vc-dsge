function mcmcredux(obj,varargin)

% mcmcredux
% 
% Combine MCMC sample chains and save part of it. saved file contains a 
% structure named "draws" with fields:
%   draw.N (number of draws)
%   draw.Param (compact matrix with draws for Param)
%   draw.LPDF (row matrix with draws for posterior log-pdf)
%   draw.AuxParam (if op.Draws.AuxParam=1, matrix with draws for AuxParam)
%
% see also:
% DSGE.Model
%
% ............................................................................
%
% Created: April 3, 2017
% Copyright (C) 2017-2018 Vasco Curdia

    %% Options
    op.Draws.BurnIn = 0.25;
    op.Draws.AuxParam = 0;
    op.Draws.CombineChains = 1;
    op.Draws.ExpandParam = 0;
    op.NDraws = 10000;
    op = updateoptions(op,varargin{:});

    fprintf('Generating MCMC Draws Redux\n')

    %% load the mcmc draws
    draws = obj.loadmcmcdraws(op.Draws,'silent',1);
    nThinning = draws.N/op.NDraws;
    idxdraws = ceil(nThinning:nThinning:draws.N);
    draws.Param = draws.Param(:,idxdraws);
    draws.LPDF = draws.LPDF(idxdraws);
    draws.N = length(idxdraws);
    fprintf('Number of draws kept: %.0f\n',draws.N);

    %% Save MCMC draws Redux
    fn = sprintf('%s-mcmc-redux',obj.Name);
    obj.Post.MCMCSample.FileNameRedux = fn;
    save(fn,'-struct','draws')
    fprintf('Saved MCMC draws redux to: %s.mat\n',fn)

end

