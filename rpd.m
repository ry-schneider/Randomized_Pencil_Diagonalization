function [S,T,D,splits,fails,flop_count] = rpd(A,B,error)
% -----------------------------------------------------------------------
% RPD wraps around DNC_EIG to compute a diagoanlization of (A,B). With 
% high probability:
%       norm(A-S*D/T,2) <= error
%       norm(B - S/T,2) <= error
%
% Requires: norm(A,2) <= 1 and norm(B,2) <= 1.
%
% D and T contain the eigenvalues and eigenvectors of (A,B) respectively.
% splits, fails, and flop_counts record performance data of the
% divide-and-conquer process, as specificed in DNC_EIG.
% 
% Calls: DNC_EIG
% -----------------------------------------------------------------------
n = size(A,1);
% --------------------------------------------------
% Set relaxed parameters (original are in comments)
% --------------------------------------------------
gamma = error/16;
alpha = 3/2; % (1/2)*ceil(2*(log(1/gamma)/log(n))+3);
epsilon = gamma/n; % gamma.^5/(64*n.^((11*alpha+25)/3)+gamma.^5);
beta = epsilon; % error*gamma.^2/(12*(1+4*gamma)*n.^(3*alpha+5)+error*gamma*n.^((2*alpha+5)/2));
omega = epsilon; % gamma.^4/(4*n.^((8*alpha+13)/3));
% ------------------
% Perturb and scale
% ------------------
G1 = (1/sqrt(2*n))*(randn(n)+1i*randn(n));
G2 = (1/sqrt(2*n))*(randn(n)+1i*randn(n));
Atilde = A+gamma*G1;
Btilde = B+gamma*G2;
% ----------------------------
% Construct the random grid g
% ----------------------------
z = complex(-4+rand*omega,-4+rand*omega);
v_left = real(z);
v_right = v_left + ceil(8/omega)*omega;
h_bottom = imag(z);
h_top = imag(z) + ceil(8/omega)*omega;
g = [v_left; v_right; h_bottom; h_top];
% --------------------------------
% Call divide-and-conquer routine
% --------------------------------
splits = [];
fails = [];
[T,D1,D2,splits,fails,flop_count] = dnc_eig(n,Atilde,n.^alpha*Btilde,epsilon,alpha,g,beta,omega,1/n,splits,fails,0);
D = n.^alpha*D1/D2;
S = Btilde*T;
end

