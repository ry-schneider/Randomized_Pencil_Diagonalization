function [P_p,P_m,g_p,g_m,n_p,n_m,fails,splits,flops] = split(A,epsilon,g,omega,beta,n,fails,splits,flops)
%--------------------------------------------------------------------------
% SPLIT subroutine defined by Banks et al. (Foundations of Computational
% Math 2022).
%
% Inputs:
%   A - m x m matrix
%   epsilon - pseudospectral guarantee
%   g - shattering grid
%   omega - box size of g
%   beta - desired projector accuracy
%   n - original problem size
%
% Additional inputs track split sizes/failed runs.
%
% Requires: epsilon pseudospectrum of A shattered w.r.t. g
% 
% Outputs:
%  P_p/P_m - spectral projectors
%  g_p/g_m - half grids corresponding to P_p/P_m
%  n_p/n_m - eigenvalue counts corresponding to g_p/g_m
%
% Calls: SGN
% -------------------------------------------------------------------------
m = size(A,1);
s_v = ceil((g(2)-g(1))/omega)+1; % number of vertical grid lines
s_h = ceil((g(4)-g(3))/omega)+1; % number of horizontal grid lines
d_g = omega*sqrt(s_v.^2 + s_h.^2);
%---------------------------------------
% Search over vertical grid lines first
%---------------------------------------
no_vertical_split = false;
line_found = false;
split_failed = false;
lines_checked = 0;
index = max(ceil(s_v/2),1);
left_index = 1;
right_index = s_v;
if index == 1
    no_vertical_split = true;
end
while line_found == false && no_vertical_split == false
    h = g(1)+(index-1)*omega;
    t = round(real(trace(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n))),0);
    lines_checked = lines_checked + 1;
    if abs(t) < min([3*m/5 m-1])
        line_found = true;
    elseif t < 0 % search to the left
        right_index = index;
        index = max(floor((left_index + index)/2),1);
    else % search to the right
        left_index = index;
        index = min(floor((right_index + index)/2),s_v); 
    end
    if right_index-left_index <= 1 || lines_checked >= floor(log2(s_v)+1)
        no_vertical_split = true;
        vertical_lines_checked = lines_checked;
    end
end
%--------------------------------------
% Apply vertical split if one is found
%--------------------------------------
if line_found == true
    disp('vertical split')
    P_p = (1/2)*(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n)+eye(m));
    P_m = (1/2)*(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n)-eye(m));
    g_p = [h;g(2);g(3);g(4)];
    g_m = [g(1);h;g(3);g(4)];
    n_p = round((t+m)/2,0);
    n_m = m-n_p;
    if m > 3
        splits = [splits; n_p/m];
    end
    flops = flops + lines_checked*m.^3;
%-------------------------------------------
% Search over horizontal grid lines instead
%-------------------------------------------
else
    lines_checked = 0;
    index = max(ceil(s_h/2),1);
    bottom_index = 1;
    top_index = s_h;
    A = -1i*A;
    while line_found == false && split_failed == false
        h = g(3)+(index-1)*omega;
        t = round(real(trace(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n))),0);
        lines_checked = lines_checked + 1;
        if abs(t) < min([3*m/5 m-1])
            line_found = true;
        elseif t < 0 % search below
            top_index = index;
            index = max(floor((index+bottom_index)/2),1);
        else % search above
            bottom_index = index;
            index = min(floor((index+top_index)/2),s_h);
        end
        if top_index-bottom_index <= 1 || lines_checked >= floor(log2(s_h)+1)
            split_failed = true;
        end
    end
    if split_failed == true
        disp('error: no dividing line found')
        fails = [fails; m];
        % ---------
        % Call EIG
        % ---------
        [V,D] = eig(-1i*A);
        P_p = V;
        P_m = D;
        g_p = g;
        g_m = g;
        n_p = m;
        n_m = m;
        flops = flops + (lines_checked+vertical_lines_checked)*m.^3;   
    %------------------------    
    % Apply horizontal split
    %------------------------
    else 
        disp('horizontal split')
        P_p = (1/2)*(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n)+eye(m));
        P_m = (1/2)*(sgn(A-h*eye(m),epsilon/4,1-(epsilon/(2*d_g.^2)),beta,n)-eye(m));
        g_p = [g(1);g(2);h;g(4)];
        g_m = [g(1);g(2);g(3);h];
        n_p = round((t+m)/2,0);
        n_m = m-n_p;
        if m > 3
            splits = [splits; n_p/m];
        end
        flops = flops + (lines_checked+vertical_lines_checked)*m.^3;
    end
end

