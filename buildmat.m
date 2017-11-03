function A = buildmat(f,x,m,n,squeezeDim2)
    
% buildmat
%
% Builds matrix A using a vector of function handles (f) and passing as 
% argument the column vector x. Matrix output is reshaped to m rows, n columns, 
% and the 3rd dimension is set to number of columns in x.
%
% Usage:
%   A = buildmat(f,x,m,n)
%
% ...........................................................................
%
% Created: March 23, 2017 by Vasco Curdia
% Copyright 2017 by Vasco Curdia

%% Option
squeezeDim2 = 0;
    
if length(f)~=m*n
    error('Length of f needs to match m*n.')
end
nx = size(x,2);
A = zeros(m,n,nx);
for jx=1:nx
    Aj = zeros(m,n);
    xj = num2cell(x(:,jx));
    for j=1:m*n
        Aj(j) = f{j}(xj{:});
    end
    A(:,:,jx) = Aj;
end
if n==1 && squeezeDim2
    A = permute(A,[1,3,2]);
end
