function [U_R,U_L] = deflate(A,B,p,k)
% -----------------------------------------------------------------------
% DEFLATE computes n x k matrices U_R and U_L whose orthonormal columns 
% span right/left deflating subspaces of (A,B) of dimension k. p is the 
% number of steps of IRS to be applied. It and k are computed externally.
%
% Outputs:
%   U_R,U_L - matrices with orthonormal columns (typically nonsquare)
%
% Calls: IRS and GRURV.
% -----------------------------------------------------------------------
[A_p,B_p] = irs(A,B,p);
[U1,~,~,~] = grurv(A_p+B_p,A_p,-1,1);
[A_p,B_p] = irs(A',B',p);
[U2,~,~,~] = grurv(A_p',(A_p+B_p)',1,-1);
U_R = U1(:,1:k);
U_L = U2(:,1:k);
end

