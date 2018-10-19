function initializeprior(obj,varargin)

% initializeprior
%
% initialize and prepare prior distribution for further use
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% 
% Copyright 2017-2018 by Vasco Curdia


%% Preamble
fprintf('\nInitializing prior\n')

%% Analyze prior distributions
np = obj.Param.N;
pNames = obj.Param.Names;
obj.Prior.Dist = obj.Param.PriorDist;
obj.Prior.Mean = obj.Param.PriorMean;
obj.Prior.SD = obj.Param.PriorSD;
obj.Prior.Mode = nan(np,1);
obj.Prior.Median = nan(np,1);
obj.Prior.Prc05 = nan(np,1);
obj.Prior.Prc95 = nan(np,1);
obj.Prior.DistParam = cell(np,1);
obj.Prior.PDFCmd = cell(np,1);
obj.Prior.RndCmd = cell(np,1);
pOptions = optimoptions(@fsolve);
pOptions.Display = 'off';
for j=1:np

    if strcmp(obj.Prior.Dist{j},'C')
        obj.Prior.Mode(j) = obj.Prior.Mean(j);
        obj.Prior.SD(j) = 0;
        obj.Prior.Median(j) = obj.Prior.Mean(j);
        obj.Prior.Prc05(j) = obj.Prior.Mean(j);
        obj.Prior.Prc95(j) = obj.Prior.Mean(j);
        obj.Prior.PDFCmd{j} = @(x)repmat(1,size(x));
        obj.Prior.RndCmd{j} = @(n)repmat(obj.Prior.Mean(j),1,n);
    
    elseif strcmp(obj.Prior.Dist{j},'N')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        obj.Prior.Mode(j) = pmean;
        obj.Prior.Median(j) = pmean;
        obj.Prior.Prc05(j) = norminv(0.05,pmean,psd);
        obj.Prior.Prc95(j) = norminv(0.95,pmean,psd);
        obj.Prior.DistParam{j} = [pmean,psd];
        obj.Prior.PDFCmd{j} = @(x)normpdf(x,pmean,psd);
        obj.Prior.RndCmd{j} = @(n)normrnd(pmean,psd,1,n);
    
    elseif strcmp(obj.Prior.Dist{j},'TN')
        % Assume x>=0
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = -pmean/psd;
        acdf = normcdf(a,0,1);
        aZ = 1-normcdf(a,0,1);
        alambda = normpdf(a,0,1)/aZ;
        adelta = alambda*(alambda-a);
        obj.Prior.Mean(j) = pmean + psd*alambda;
        obj.Prior.SD(j) = psd*(1-adelta)^(1/2);
        obj.Prior.Mode(j) = max(0,pmean);
        obj.Prior.Median(j) = norminv(0.5*aZ+acdf,pmean,psd);
        obj.Prior.Prc05(j) = norminv(0.05*aZ+acdf,pmean,psd);
        obj.Prior.Prc95(j) = norminv(0.95*aZ+acdf,pmean,psd);
        obj.Prior.DistParam{j} = [pmean,psd];
        obj.Prior.PDFCmd{j} = @(x)(x>=0).*normpdf((x-pmean)/psd,0,1)/psd/aZ;
        obj.Prior.RndCmd{j} = @(n)norminv(rand(1,n).*aZ+zcdf,pmean,psd);
    
    elseif strcmp(obj.Prior.Dist{j},'B')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = pmean*(pmean-pmean^2-psd^2)/psd^2;
        b = a*(1/pmean-1);
        obj.Prior.Mode(j) = min(max(0,(a-1)/(a+b-2)),1);
        obj.Prior.Median(j) = betainv(0.5,a,b);
        obj.Prior.Prc05(j) = betainv(0.05,a,b);
        obj.Prior.Prc95(j) = betainv(0.95,a,b);
        obj.Prior.DistParam{j} = [a,b];
        obj.Prior.PDFCmd{j} = @(x)betapdf(x,a,b);
        obj.Prior.RndCmd{j} = @(n)betarnd(a,b,1,n);
        
    elseif strcmp(obj.Prior.Dist{j},'G')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        a = (pmean/psd)^2;
        b = pmean/a;
        if a>=1
            obj.Prior.Mode(j) = (a-1)*b;
        else
            obj.Prior.Mode(j) = NaN;
        end
        obj.Prior.Median(j) = gaminv(0.5,a,b);
        obj.Prior.Prc05(j) = gaminv(0.05,a,b);
        obj.Prior.Prc95(j) = gaminv(0.95,a,b);
        obj.Prior.DistParam{j} = [a,b];
        obj.Prior.PDFCmd{j} = @(x)gampdf(x,a,b);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b,1,n);
        
    elseif strcmp(obj.Prior.Dist{j},'IG1')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        if psd==inf
            a = 1;
        else
%             [a,rc] = csolvevb(@(x)igamsolve(x,pmean,psd),5,[],1e-10,1000);
            [a,~,rc] = fsolve(@(x)igamsolve(x,pmean,psd),5,pOptions);
%             if rc~=0, 
            if rc~=1, 
                error('Search for iGam parameters failed, rc=%.0f',rc), 
            end
        end
        b = (gamma(a-1/2)/pmean/gamma(a))^2;
        obj.Prior.Mode(j) = (1/b/(a+1/2))^(1/2);
        obj.Prior.Median(j) = gaminv(0.5,a,b)^(-1/2);
        obj.Prior.Prc05(j) = gaminv(1-0.05,a,b)^(-1/2);
        obj.Prior.Prc95(j) = gaminv(1-0.95,a,b)^(-1/2);
        obj.Prior.DistParam{j} = [a,b];
        obj.Prior.PDFCmd{j} = @(x)(x>0).*(gampdf(x.^(-2),a,b)*2./x.^3);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b,1,n).^(-1/2);
    
    elseif strcmp(obj.Prior.Dist{j},'IG2')
        pmean = obj.Prior.Mean(j);
        psd = obj.Prior.SD(j);
        if psd==inf
            a = 2;
        else
            a = 2+pmean^2/psd^2;
        end
        b = 1/pmean/(a-1);
        obj.Prior.Mode(j) = 1/b/(a+1);
        obj.Prior.Median(j) = gaminv(0.5,a,b)^(-1);
        obj.Prior.Prc05(j) = gaminv(1-0.05,a,b)^(-1);
        obj.Prior.Prc95(j) = gaminv(1-0.95,a,b)^(-1);
        obj.Prior.DistParam{j} = [a,b];
        obj.Prior.PDFCmd{j} = @(x)(x>0).*(gampdf(x.^(-1),a,b)./x.^2);
        obj.Prior.RndCmd{j} = @(n)gamrnd(a,b).^(-1);
        
    end
end

obj.analyzeprior(varargin{:})

end

%% functions used

function f = igamsolve(a,pmean,psd)
    f = 1./(a-1).*(pmean*gamma(a)./gamma(a-1/2)).^2-pmean^2-psd^2;
end
