function analyzedist(obj)

% analyzedist
%
% Analyzes prior distribution
%
% See also:
% DSGE.Prior
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% 
% Copyright 2017 by Vasco Curdia


%% Preamble

%% Analyze Parameters
np = obj.NParam;
pNames = obj.ParamNames;
obj.Mode = nan(np,1);
obj.Median = nan(np,1);
obj.Prc05 = nan(np,1);
obj.Prc95 = nan(np,1);
obj.DistParam = cell(np,1);
obj.LPdfCmd = cell(np,1);
obj.PdfCmd = cell(np,1);
obj.RndCmd = cell(np,1);
pOptions = optimoptions(@fsolve);
pOptions.Display = 'off';
for j=1:np

    if strcmp(obj.Dist{j},'C')
        obj.Mode(j) = obj.Mean(j);
        obj.SD(j) = 0;
        obj.Median(j) = obj.Mean(j);
        obj.Prc05(j) = obj.Mean(j);
        obj.Prc95(j) = obj.Mean(j);
        obj.LPdfCmd{j} = 0;
        obj.PdfCmd{j} = 1;
        obj.RndCmd{j} = @(n)repmat(obj.Mean(j),1,n);
    
    elseif strcmp(obj.Dist{j},'N')
        pmean = obj.Mean(j);
        psd = obj.SD(j);
        obj.Mode(j) = pmean;
        obj.Median(j) = pmean;
        obj.Prc05(j) = norminv(0.05,pmean,psd);
        obj.Prc95(j) = norminv(0.95,pmean,psd);
        obj.DistParam{j} = [pmean,psd];
        obj.LPdfCmd{j} = sprintf('log(normpdf(%s,%.16f,%.16f))',...
                                    pNames{j},pmean,psd);
        obj.PdfCmd{j} = sprintf('normpdf(%%.16f,%.16f,%.16f)',pmean,psd);
        obj.RndCmd{j} = @(n)normrnd(pmean,psd,1,n);
    
    elseif strcmp(obj.Dist{j},'TN')
        % Assume x>=0
        pmean = obj.Mean(j);
        psd = obj.SD(j);
        a = -pmean/psd;
        acdf = normcdf(a,0,1);
        aZ = 1-normcdf(a,0,1);
        alambda = normpdf(a,0,1)/aZ;
        adelta = alambda*(alambda-a);
        obj.Mean(j) = pmean + psd*alambda;
        obj.SD(j) = psd*(1-adelta)^(1/2);
        obj.Mode(j) = max(0,pmean);
        obj.Median(j) = norminv(0.5*aZ+acdf,pmean,psd);
        obj.Prc05(j) = norminv(0.05*aZ+acdf,pmean,psd);
        obj.Prc95(j) = norminv(0.95*aZ+acdf,pmean,psd);
        obj.DistParam{j} = [pmean,psd];
        obj.LPdfCmd{j} = sprintf(...
            ['log((%1$s>=0)'...
             '*normpdf((%1$s-%2$.16f)/%3$.16f,0,1)/%3$.16f/%4$.16f)'],...
            pNames{j},pmean,psd,aZ);
        obj.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>=0)',...
             '*normpdf((%%1$.16f-%1$.16f)/%2$.16f,0,1)/%2$.16f/%3$.16f'],...
            pmean,psd,aZ);
        obj.RndCmd{j} = @(n)norminv(rand(1,n).*aZ+zcdf,pmean,psd);
    
    elseif strcmp(obj.Dist{j},'B')
        pmean = obj.Mean(j);
        psd = obj.SD(j);
        a = pmean*(pmean-pmean^2-psd^2)/psd^2;
        b = a*(1/pmean-1);
        obj.Mode(j) = min(max(0,(a-1)/(a+b-2)),1);
        obj.Median(j) = betainv(0.5,a,b);
        obj.Prc05(j) = betainv(0.05,a,b);
        obj.Prc95(j) = betainv(0.95,a,b);
        obj.DistParam{j} = [a,b];
        obj.LPdfCmd{j} = sprintf('log(betapdf(%s,%.16f,%.16f))',pNames{j},a,b);
        obj.PdfCmd{j} = sprintf('betapdf(%%.16f,%.16f,%.16f)',a,b);
        obj.RndCmd{j} = @(n)betarnd(a,b,1,n);
        
    elseif strcmp(obj.Dist{j},'G')
        pmean = obj.Mean(j);
        psd = obj.SD(j);
        a = (pmean/psd)^2;
        b = pmean/a;
        if a>=1
            obj.Mode(j) = (a-1)*b;
        else
            obj.Mode(j) = NaN;
        end
        obj.Median(j) = gaminv(0.5,a,b);
        obj.Prc05(j) = gaminv(0.05,a,b);
        obj.Prc95(j) = gaminv(0.95,a,b);
        obj.DistParam{j} = [a,b];
        obj.LPdfCmd{j} = sprintf('log(gampdf(%s,%.16f,%.16f))',pNames{j},a,b);
        obj.PdfCmd{j} = sprintf('gampdf(%%.16f,%.16f,%.16f)',a,b);
        obj.RndCmd{j} = @(n)gamrnd(a,b,1,n);
        
    elseif strcmp(obj.Dist{j},'IG1')
        pmean = obj.Mean(j);
        psd = obj.SD(j);
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
        obj.Mode(j) = (1/b/(a+1/2))^(1/2);
        obj.Median(j) = gaminv(0.5,a,b)^(-1/2);
        obj.Prc05(j) = gaminv(0.05,a,b)^(-1/2);
        obj.Prc95(j) = gaminv(0.95,a,b)^(-1/2);
        obj.DistParam{j} = [a,b];
        obj.LPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-2),%2$.16f,%3$.16f)*2/%1$s^3))',...
            pNames{j},a,b);
        obj.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-2),%1$.16f,%2$.16f)*2/%%1$.16f^3)'],a,b);
        obj.RndCmd{j} = @(n)gamrnd(a,b,1,n).^(-1/2);
    
    elseif strcmp(obj.Dist{j},'IG2')
        pmean = obj.Mean(j);
        psd = obj.SD(j);
        if psd==inf
            a = 2;
        else
            a = 2+pmean^2/psd^2;
        end
        b = 1/pmean/(a-1);
        obj.Mode(j) = 1/b/(a+1);
        obj.Median(j) = gaminv(0.5,a,b)^(-1);
        obj.Prc05(j) = gaminv(0.05,a,b)^(-1);
        obj.Prc95(j) = gaminv(0.95,a,b)^(-1);
        obj.DistParam{j} = [a,b];
        obj.LPdfCmd{j} = sprintf(...
            'log((%1$s>0)*(gampdf(%1$s^(-1),%2$.16f,%3$.16f)/%1$s^2))',...
            pNames{j},a,b);
        obj.PdfCmd{j} = sprintf(...
            ['(%%1$.16f>0)',...
             '*(gampdf(%%1$.16f^(-1),%1$.16f,%2$.16f)/%%1$.16f^2)'],a,b);
        obj.RndCmd{j} = @(n)gamrnd(a,b).^(-1);
        
    end
end


end

%% functions used

function f = igamsolve(a,pmean,psd)
    f = 1./(a-1).*(pmean*gamma(a)./gamma(a-1/2)).^2-pmean^2-psd^2;
end
