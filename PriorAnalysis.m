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
    PrcList{jprc} = sprintf('PriorPrc%03.0f',1000*op.Percentiles(jprc));
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
    elseif strcmp(p.PriorDist{j},'TN')
        % Assume x>=0
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        a = -pmean/pse;
        acdf = normcdf(a,0,1);
        aZ = 1-normcdf(a,0,1);
        alambda = normpdf(a,0,1)/aZ;
        adelta = alambda*(alambda-a);
        p.PriorMean(j) = pmean + pse*alambda;
        p.PriorSE(j) = pse*(1-adelta)^(1/2);
        p.PriorMode(j) = max(0,pmean);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = norminv(...
                op.Percentiles(jprc)*aZ+acdf,pmean,pse);
        end
        p.PriorParams(j,:) = [pmean,pse];
        p.PriorLPdfCmd{j} = sprintf(...
            ['log((%1$s>=0)'...
             '*normpdf((%1$s-%2$.16f)/%3$.16f,0,1)/%3$.16f/%4$.16f)'],...
            p.Names{j},pmean,pse,aZ);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>=0)',...
             '*normpdf((%%1$.16f-%1$.16f)/%2$.16f,0,1)/%2$.16f/%3$.16f'],...
            pmean,pse,aZ);
        p.PriorRndCmd{j} = sprintf('norminv(rand*%.16f+%.16f,%.16f,%.16f)',...
                                   aZ,acdf,pmean,pse);
    elseif strcmp(p.PriorDist{j},'B')
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        a = pmean*(pmean-pmean^2-pse^2)/pse^2;
        b = a*(1/pmean-1);
        p.PriorMode(j) = min(max(0,(a-1)/(a+b-2)),1);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = betainv(op.Percentiles(jprc),a,b);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf('log(betapdf(%s,%.16f,%.16f))',...
                                    p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf('betapdf(%%.16f,%.16f,%.16f)',a,b);
        p.PriorRndCmd{j} = sprintf('betarnd(%.16f,%.16f)',a,b);
    elseif strcmp(p.PriorDist{j},'G')
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        a = (pmean/pse)^2;
        b = pmean/a;
        if a>=1
            p.PriorMode(j) = (a-1)*b;
        else
            p.PriorMode(j) = NaN;
        end
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(op.Percentiles(jprc),a,b);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf('log(gampdf(%s,%.16f,%.16f))',...
                                    p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf('gampdf(%%.16f,%.16f,%.16f)',a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)',a,b);
    elseif strcmp(p.PriorDist{j},'IG1')
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        if pse==inf
            a = 1;
        else
            fname = sprintf('igamsolve%.0f',cputime*1e10);
            fid=fopen([fname,'.m'],'wt');
            fprintf(fid,'function f=%s(x)\n',fname);
            fprintf(fid,'pmean = %.16f;\n',pmean);
            fprintf(fid,'pvar = %.16f;\n',pse^2);
            fprintf(fid,'for j=1:length(x)\n');
            fprintf(fid,'    a = x(j);\n');
            fprintf(fid,...
                    ['    f(j) = 1/(a-1)*(pmean*gamma(a)/gamma(a-1/2))^2',...
                     '-pmean^2-pvar;\n']);
            fprintf(fid,'end\n');
            fclose(fid);
            [a,rc] = csolvevb(fname,5,[],1e-10,1000);
            if rc~=0, 
                error('Search for iGam parameters failed, rc=%.0f',rc), 
            end
            delete([fname,'.m'])
        end
        b = (gamma(a-1/2)/pmean/gamma(a))^2;
        p.PriorMode(j) = (1/b/(a+1/2))^(1/2);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(1-op.Percentiles(jprc),a,b)^(-1/2);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-2),%2$.16f,%3$.16f)*2/%1$s^3))',...
            p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-2),%1$.16f,%2$.16f)*2/%%1$.16f^3)'],a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)^(-1/2)',a,b);
    elseif strcmp(p.PriorDist{j},'IG2')
        pmean = p.PriorMean(j);
        pse = p.PriorSE(j);
        if pse==inf
            a = 2;
        else
            a = 2+pmean^2/pse^2;
        end
        b = 1/pmean/(a-1);
        p.PriorMode(j) = 1/b/(a+1);
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = gaminv(1-op.Percentiles(jprc),a,b)^(-1);
        end
        p.PriorParams(j,:) = [a,b];
        p.PriorLPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',...
            p.Names{j},a,b);
        p.PriorPdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)'],a,b);
        p.PriorRndCmd{j} = sprintf('gamrnd(%.16f,%.16f)^(-1)',a,b);
    elseif strcmp(p.PriorDist{j},'C')
        pmean = p.PriorMean(j);
        pse = 0;
        p.PriorSE(j) = 0;
        p.PriorMode(j) = pmean;
        for jprc=1:nprc
            p.(PrcList{jprc})(j) = pmean;
        end
        p.PriorParams(j,:) = [pmean,pse];
        p.PriorLPdfCmd{j} = sprintf('log(1*(%s==%.16f))',p.Names{j},pmean);
        p.PriorPdfCmd{j} = sprintf('1*(%%.16f==%.16f)',pmean);
        p.PriorRndCmd{j} = sprintf('%.16f',pmean);
    end
end
s.Param = p;


%% display results on screen
if op.ShowTable
    fprintf('\nPrior Information:')
    fprintf('\n------------------\n')
    namelength = [cellfun('length',p.Names)];
    namelengthmax = max(namelength);
    DispList = {'','','Names';
                'Prior','dist','PriorDist';
                '','  mode','PriorMode';
                '','  mean','PriorMean';
                '','   se','PriorSE';
                '','   5%','PriorPrc050';
                '',' median','PriorPrc500';
                '','   95%','PriorPrc950';
               }';
    nc = size(DispList,2);
    for jr=1:2
        fprintf(['%-',int2str(namelengthmax),'s'],DispList{jr,1});
        fprintf('  %-4s',DispList{jr,2});
        for jc=3:nc
            fprintf('  %-7s',DispList{jr,jc});
        end
        fprintf('\n');
    end
    for j=1:np
        fprintf(['%',int2str(namelengthmax),'s'],p.(DispList{3,1}){j});
        fprintf('  %4s',p.(DispList{3,2}){j});
        for jc=3:nc
            fprintf('  %7.4f',p.(DispList{3,jc})(j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end


%% Generate prior draws
fprintf('Generating function to draw from prior...\n')
s.FileName.PriorDraw = [s.Spec,'PriorDraw'];


%% -------------------------------------------------------------------

%% Finish up
s.Status.(Action) = 1;
s.TimeElapsed.(Action) = toc-s.TimeElapsed.(Action);

%% -------------------------------------------------------------------

