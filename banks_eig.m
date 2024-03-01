function [V,D,fails,splits,flops] = banks_eig(A,delta,g,omega,epsilon,theta,n,fails,splits,flops,stop)
%--------------------------------------------------------------------------
% EIG subroutine defined by Banks et al. (Foundations of Computational 
% Math 2022)
%
% Inputs:
%   A - m x m matrix
%   delta - eigenvector accuracy
%   g - shattering grid
%   omega - box size of g
%   epsilon - pseudospectral guarantee
%   theta - failure probability
%   n - original problem size
% 
% The input stop determines the problem size at which divide-and-conquer
% passess off to QR. Other inputs track split size/failed runs.
%
% Requires:
%   1. epsilon pseudospectrum of A shattered with respect to g
%   2. m <= n
%
% Outputs:
%   V - eigenvector matrix
%   D - diagaonal eigenvalue matrix
%
% Calls: SPLIT, BANKS_DEFLATE, and BANKS_EIG
%--------------------------------------------------------------------------
m = size(A,1);
%-------------------
% Stopping criteria
%-------------------
if m == 1 
    V = 1;
    D = A;
elseif m <= stop
    [V,D] = eig(A);
    for i = 1:m
        V(:,i) = V(:,i)/norm(V(:,i)/2);
    end
%------------------------------------
% Run one step of divide-and-conquer
%------------------------------------
else
    eta = delta*epsilon.^2/200;
    beta = eta.^4*theta.^2/((20*n).^6*4*n.^8);
    fail_count = size(fails,1);
    [P_p,P_m,g_p,g_m,n_p,n_m,fails,splits,flops] = split(A,epsilon,g,omega,beta,n,fails,splits,flops);
    if size(fails,1) > fail_count
        V = P_p;
        D = P_m;
    elseif n_p+n_m ~= m
        disp('error: split computed inccorrectly')
        [V,D] = eig(A);
        fails = [fails; 0];
    else
        Q_p = banks_deflate(P_p,n_p);
        Q_m = banks_deflate(P_m,n_m);
        A_p = Q_p'*A*Q_p;
        A_m = Q_m'*A*Q_m;
        [V_p,D_p,fails,splits,flops] = banks_eig(A_p,4*delta/5,g_p,omega,4*epsilon/5,theta,n,fails,splits,flops,stop);
        [V_m,D_m,fails,splits,flops] = banks_eig(A_m,4*delta/5,g_m,omega,4*epsilon/5,theta,n,fails,splits,flops,stop);
        D = [D_p zeros(n_p,n_m); zeros(n_m,n_p) D_m];
        V = [Q_p*V_p Q_m*V_m];
    end
    V = V/norm(V,2);
end
end

