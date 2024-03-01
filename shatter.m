function [X,g,omega,epsilon] = shatter(A,gamma)
%--------------------------------------------------------------------------
% SHATTER subroutine defined by Banks et al. (Foundations of 
% Computational Math 2022)
%
% Inputs: A (with norm at most one) and 0 < gamma < 1/2
% Outputs:
%   X - perturbed matrix
%   g - shattering grid (stored as a vector of four real numbers)
%   omega - box size of g
%   epsilon - shattering epsilon
%--------------------------------------------------------------------------
n = size(A,1);
% ---------------------
% Draw G and perturb A
% ---------------------
G = (1/sqrt(2*n))*(randn(n)+1i*randn(n));
X = A+gamma*G;
% ----------------------------
% Construct the random grid g
% ----------------------------
omega = gamma/n; % gamma.^4/(4*n.^5);
z = complex(-4+rand*omega,-4+rand*omega);
v_left = real(z); % left most vertical grid line
v_right = v_left + ceil(8/omega)*omega; % right most vertical grid line
h_bottom = imag(z); % bottom most horizontal grid line
h_top = imag(z) + ceil(8/omega)*omega; % top most horizontal grid line
g = [v_left; v_right; h_bottom; h_top];
% ----------------------
% Set shattered epsilon
% ----------------------
epsilon = gamma/n; % gamma.^4/(16*n.^9);
end

