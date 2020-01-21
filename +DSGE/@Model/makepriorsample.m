function makepriorsample(obj,varargin)

% makepriorsample
%
% makes a sample from the prior distribution
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: October 19, 2018 by Vasco Curdia
% Copyright 2018 by Vasco Curdia


%% Preamble
fprintf('Making prior sample\n')

%% Options
op.NDraws = 10000;
op.Percentiles = [0.01, 0.05, 0.15, 0.25, 0.75, 0.85, 0.95, 0.99];
op = updateoptions(op,varargin{:});

%% Make draws
BadDraws = false(1,op.NDraws);
xd = nan(obj.Param.N,op.NDraws);
for j=1:obj.Param.N
    xd(j,:) = obj.Prior.RndCmd{j}(op.NDraws);
end
xdAux = zeros(obj.AuxParam.N,op.NDraws);
parfor jd=1:op.NDraws
    Matsj = obj.solveree(xd(:,jd));
    BadDraws(jd) = ~Matsj.Status;
    xdAux(:,jd) = Matsj.AuxParam;
end
obj.Prior.Sample.NDraws = op.NDraws;
obj.Prior.Sample.NBadDraws = sum(BadDraws);
obj.Prior.Sample.FractionBadDraws = obj.Prior.Sample.NBadDraws/op.NDraws;
obj.Prior.LPDFCorrection = -log(1-obj.Prior.Sample.FractionBadDraws);
xd(:,BadDraws) = [];
xdAux(:,BadDraws) = [];
obj.Prior.Sample.NDrawsUsed = size(xd,2);
fprintf('Number of accepted draws: %.0f\n',obj.Prior.Sample.NDrawsUsed);
fprintf('Percent of rejected draws: %.2f%%\n',...
        obj.Prior.Sample.FractionBadDraws*100);
fprintf('log-prior correction: %.6f\n',obj.Prior.LPDFCorrection);

obj.Prior.Sample.Param = sumstats(xd,op.Percentiles);
obj.Prior.Sample.AuxParam = sumstats(xdAux,op.Percentiles);

% Save prior sample
fn = sprintf('%s-priorsample',obj.Name);
obj.Prior.Sample.FileName = fn;
draws.N = obj.Prior.Sample.NDrawsUsed;
draws.Param = xd;
save(fn,'-struct','draws')
fprintf('Saved prior sample to: %s.mat\n',fn)


