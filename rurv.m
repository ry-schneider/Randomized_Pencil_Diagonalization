function [U,R,V] = rurv(A)
% --------------------------------------------------------------------
% RURV computes a randomized, rank-revealing factorization A = U*R*V.
%
% Outputs:
%   U - unitary
%   R - upper triangular
%   V - Haar
%
% See Demmel, Dumitriu and Holtz (Numerische Mathematik 2007) for
% rank-revealing guarantees.
% --------------------------------------------------------------------
n = size(A,1);
B = (1/sqrt(2))*(randn(n)+1i*randn(n));
[V,~] = qr(B);
[U,R] = qr(A*V');
end

