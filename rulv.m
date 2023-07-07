function [U,L,V] = rulv(A)
% -----------------------------------------------------------------------
% RULV computes a randomized rank-revealing factorization A = U*L*V
%  
% Outputs:
%   U - unitary
%   L - lower triangular
%   V - Haar
%
% See Demmel, Dumitriu, and Holtz (Numerische Mathematik 2007).
% -----------------------------------------------------------------------
n = size(A,1);
B = (1/sqrt(2))*(randn(n)+1i*randn(n));
[V,~] = qr(B);
[R,Q] = rq(V*A');
U = Q';
L = R';
end

