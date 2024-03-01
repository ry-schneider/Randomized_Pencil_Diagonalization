function [V,D,fails, splits,flops] = banks_algorithm(A,error)
%--------------------------------------------------------------------------
% Randomized diagonalization algorithm defined by Banks et al. (Foundations
% of Computational Math 2022).
%
% Inputs:
%   A - n x n matrix
%   error - desired backward diagonalization error
% 
% Requires: norm(A,2) <= 1
% 
% Output: Backward diagonalization A = V*D*V^{-1} and performance data
% 
% Calls: SHATTER, BANKS_EIG
%--------------------------------------------------------------------------
n = size(A,1);
[X,g,omega,epsilon] = shatter(A,error/8);
fails = [];
splits = [];
delta = epsilon; % error.^3/(1536*n.^(9/2));
[V,D,fails,splits,flops] = banks_eig(X,delta,g,omega,epsilon,1/n,n,fails,splits,0);
end

