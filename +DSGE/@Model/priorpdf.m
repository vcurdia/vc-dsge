function p = priorpdf(obj,x)

% priorpdf
%
% evaluate the prior pdf for a set of parameters 
%
% See also:
% DSGE.Model
%
% ...........................................................................
% 
% Created: March 19, 2017 by Vasco Curdia
% Copyright 2017-2018 by Vasco Curdia


p = repmat(exp(obj.Prior.LPDFCorrection),1,size(x,2));
for j=1:obj.Param.N
    p = p.*obj.Prior.PDFCmd{j}(x(j,:));
end
