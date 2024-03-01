function [Q] = banks_deflate(P,k)
%--------------------------------------------------------------------------
% DEFLATE subroutine defined by Banks et al. (Foundations of
% Computational Math 2022).
%
% Inputs: 
%   P - approximate projector
%   k - approximate rank of P
%
% Outputs:
%   Q - n x k matrix containing a basis of range(P)
%
% Calls: RURV
%--------------------------------------------------------------------------
[U,~,~] = rurv(P);
Q = U(:,1:k);
end

