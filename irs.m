function [A_p,B_p] = irs(A,B,p)
% ------------------------------------------------------------------------
% IRS performs p steps of implicit repeated squaring to (A,B).
%
%   inv(B_p)*A_p = (inv(B)*A)^(2.^p)
%
% See work of Bai, Demmel, and Gu (Numerische Mathematik 1997).
% ------------------------------------------------------------------------
n = size(A,1);
A_p = A;
B_p = B;
for i = 1:p
    [Q,~] = qr([B_p; -A_p]);
    A_p = Q(1:n,(n+1):(2*n))'*A_p;
    B_p = Q((n+1):(2*n),(n+1):(2*n))'*B_p;
end
end

