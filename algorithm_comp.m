function [e_1,d_1,fails1,splits1,flops1,cond_1,qz_error,e_2,d_2,fails2,splits2,flops2,cond_2,qr_error] = algorithm_comp(A,B,error,stop)
%--------------------------------------------------------------------------
% On the input pencil (A,B), this function compares RPD with the following 
% three-step algorithm:
%   1. Perturb the input matrices A,B to obatin Atilde, Btilde
%   2. Form the product Btilde\Atilde
%   3. Call BANKS_EIG 
%
% Inputs:
%   A,B - n x n pencils
%   error - desired backward diagonalization error
%   stop - divide-and-conquer stopping criteria
% 
% Calls: DNC_EIG and BANKS_EIG
% -------------------------------------------------------------------------
n = size(A,1);
true_eigs = eig(A,B);
c = max(norm(A,2),norm(B,2));
A = A/c;
B = B/c;
%------------------------------------------------------------------
% Construct the random shattering grid g and set shattered epsilon
%------------------------------------------------------------------
gamma = error/16;
omega = gamma/n; % gamma.^4/(4*n.^5);
z = complex(-4+rand*omega,-4+rand*omega);
v_left = real(z); % left most vertical grid line
v_right = v_left + ceil(8/omega)*omega; % right most vertical grid line
h_bottom = imag(z); % bottom most horizontal grid line
h_top = imag(z) + ceil(8/omega)*omega; % top most horizontal grid line
g = [v_left; v_right; h_bottom; h_top];
epsilon = omega; %gamma.^5/(64*n.^((11*alpha+25)/3)+gamma.^5); % gamma.^5/(16*n.^((11*alpha+25)/3)); % or gamma.^5/(64*n.^((11*alpha+25)/3)+gamma.^5)
%----------------------
% Perturb the matrices 
%----------------------
G1 = (1/sqrt(2*n))*(randn(n)+1i*randn(n));
G2 = (1/sqrt(2*n))*(randn(n)+1i*randn(n));
Atilde = A+gamma*G1;
Btilde = B+gamma*G2;
%------------------------
% Apply inverse-free RPD
%------------------------
splits1 = [];
fails1 = [];
qz_eigs = eig(Atilde,Btilde);
[T,D1,D2,splits1,fails1,flops1] = dnc_eig(n,Atilde,Btilde,epsilon,0,g,epsilon,omega,1/n,splits1,fails1,0,stop);
D = zeros(n);
for i = 1:n
    D(i,i) = D1(i,i)/D2(i,i);
end
S = Btilde*T;
d_1 = max(norm(A-S*D/T,2),norm(B-S/T,2));
eig_if = diag(D);
cond_1 = cond(T,2);
%--------------------------------------
% Form the product and apply BANKS_EIG
%--------------------------------------
splits2 = [];
fails2 = [];
X = Btilde\Atilde;
qr_eigs = eig(X);
[V_banks,D_banks,fails2,splits2,flops2] = banks_eig(X,epsilon,g,omega,epsilon,1/n,n,fails2,splits2,0,stop);
S_banks = Btilde*V_banks;
d_2 = max(norm(A-S_banks*D_banks/V_banks,2),norm(B-S_banks/V_banks,2));
eig_banks = diag(D_banks);
cond_2 = cond(V_banks,2);
%---------------------------
% Compute eigenvalue errors
%---------------------------
true_eigs = sort(true_eigs,'ComparisonMethod','abs');
eig_if = sort(eig_if,'ComparisonMethod','abs');
eig_banks = sort(eig_banks,'ComparisonMethod','abs');
qz_eigs = sort(qz_eigs,'ComparisonMethod','abs');
qr_eigs = sort(qr_eigs,'ComparisonMethod','abs');
e_1 = mean(abs(true_eigs(1:n-1)-eig_if(1:n-1)));
e_2 = mean(abs(true_eigs(1:n-1)-eig_banks(1:n-1)));
qz_error = mean(abs(true_eigs(1:n-1)-qz_eigs(1:n-1)));
qr_error = mean(abs(true_eigs(1:n-1)-qr_eigs(1:n-1)));
