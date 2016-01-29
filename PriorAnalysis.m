function s = PriorAnalysis(s)

% PriorAnalysis
%
% Analyzes the priors
%
% See also:
%
% ...........................................................................
% 
% Created: January 28, 2016 by Vasco Curdia
% 
% Copyright 2016 by Vasco Curdia

%% -------------------------------------------------------------------

%% Preamble

Action = 'PriorAnalysis';

% % check if model already prepared
% if isfield(s.Status,Action) && s.Status.(Action), return, end

fprintf('\n*** Analyzing DSGE Prior distribution\n')

% Set Timer
s.TimeElapsed.(Action) = toc();

%% Options
if isfield(s.Options,'Prior')
    op = s.Options.Prior;
else
    op = struct; 
end
if ~isfield(op,'ShowTable'), op.ShowTable = 1; end
if ~isfield(op,'nDrawsSample'),op.nDrawsSample = 1000; end
if ~isfield(op,'ShowTableSample'), op.ShowTableSample = 1; end
if ~isfield(op,'Percentiles')
    op.Percentiles = [0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99];
end
s.Options.Prior = op;
s.FileName.PriorDraws = [s.Spec,'PriorDraws'];

%% -------------------------------------------------------------------

%% Prepare variables
nprc = length(op.Percentiles);
PrcList = cell(nprc,1);
for jprc=1:nprc
    PrcList{jprc} = sprintf('p.PriorPrc%03.0f',1000*op.Percentiles(jprc));
end


%% Analyze Parameters
p = s.Param;
np = s.n.Param;
p.PriorMode = nan(np,1);
p.PriorParams = nan(np,2);
for jprc=1:nprc
    p.(PrcList{jprc}) = nan(np,1);
end
p.PriorLPdfCmd = cell(np,1);
p.PriorPdfCmd = cell(np,1);
p.PriorRndCmd = cell(np,1);
for j=1:np
    if strcmp(p.PriorDist{j},'N')
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        p.PriorMode(j) = pmean;
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = norminv(op.Percentiles(jprc),pmean,pse);
        end
        p.PriorParams(j,:) = [pmean,pse];
        p.PriorLPdfCmd{j} = sprintf('log(normpdf(%s,%.16f,%.16f))',...
                                    p.Names{j},pmean,pse);
        p.PriorPdfCmd{j} = sprintf('normpdf(%%.16f,%.16f,%.16f)',pmean,pse);
        p.PriorRndCmd{j} = sprintf('normrnd(%.16f,%.16f)',pmean,pse);
    end
end
s.Param = p;

%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------

