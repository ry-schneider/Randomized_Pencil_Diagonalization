function [U,R_1,R_2,V] = grurv(A_1,A_2,m_1,m_2)
% -----------------------------------------------------------------------
% GRURV implicitly computes a rank-revealing factorizaton 
%           A_1^m_1*A_2^m_2 = UR_1^m_1*R_2^m_2*V
% for m_1,m_2 in {+1,-1}.
%
% Outputs:
%   U - unitary
%   R_1,R_2 - upper triangular
%   V - Haar
%
% This is a special two-matrix version of the generalization of RURV
% developed by Ballard, Demmel, and Dumitriu (Technical Report 2010).
%
% Calls: RURV and RULV.
% -----------------------------------------------------------------------
if m_2 == 1
    [U,R_2,V] = rurv(A_2);
else
    [U,L_2,V] = rulv(A_2');
    R_2 = L_2';
end 
Ucurrent = U;
if m_1 == 1
    [U,R_1] = qr(A_1*Ucurrent);
    Ucurrent = U;
else
    [R_1,U] = rq(Ucurrent'*A_1);
    Ucurrent = U';
end
U = Ucurrent;
end

