function [S] = sgn(A,epsilon,alpha,delta,n)
%--------------------------------------------------------------------------
% SGN subroutine defined by Banks et al. (Foundations of Computational
% Math 2022).
%
% Inputs:
%   A - m x m matrix
%   epsilon - pseudospectral guarantee
%   alpha - circle parameter
%   delta - desired accuracy
%   n - original problem size
%  
% Requires: epsilon pseudospectrum of A contained in C_alpha
%
% Outputs:
%    S - approximate sign function (with ||S - sgn(A)|| <= delta)
%--------------------------------------------------------------------------
N = ceil(log2(n/(4*epsilon))); %ceil(log2(1/(1-alpha)) + 3*log2(log2(1/(1-alpha))) + log2(log2(1/(delta*epsilon))) +  7.59);
A_j = A;
for i = 1:N
   A_j = (1/2)*(A_j+inv(A_j));
end
S = A_j;
end

