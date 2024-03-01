function [T,D1,D2,splits,fails,flop_count] = dnc_eig(n,A,B,epsilon,alpha,g,beta,omega,theta,splits,fails,flop_count,stop)
% ------------------------------------------------------------------------
% DNC_EIG executes divide-and-conquer on the pencil (A,n^(alpha)*B) over 
% the grid g. beta and theta are respectively an eigenvector accuracy
% parameter and an acceptable failure probability. 
%
% The routine requires that the epsilon-pseudospectrum of (A,n.^(alpha)*B) 
% is shattered with respect to the grid g, as defined by  Banks et al. 
% (Foundations of Computational Math 2022).
%
% The grid g is stored as a vector of four real numbers; the first two 
% entries record the real parts of the left/right-most vertical grid lines 
% while the second two record the imaginary parts of the bottom/top-most 
% horizontal lines. All other grid lines can be obtained from these given
% the box size omega. Note that g need not be square.
%
% Outputs:
%   T - matrix of right eigenvectors
%   (D1,D2) - diagonal pencil containing the eigenvalues of (A,n^(alpha)*B)
%   splits - array that records the split size at each step 
%   fails - array that records the problem sizes for default calls to QZ 
%           when they occur
%   flop_count - pseudo-flop count (at each step the size of the problem
%                cubed times the number of grid lines checked)
%
% NOTE: DNC_EIG CALLS ITSELF RECURSIVELY
% Additional calls: DEFLATE, GRURV, IRS
% ------------------------------------------------------------------------
m = size(A,1);
s_v = ceil((g(2)-g(1))/omega)+1;
s_h = ceil((g(4)-g(3))/omega)+1;
s = max(s_v, s_h);
% -----------------------------
% Recursion stopping condition
% -----------------------------
if m == 1 
    T = 1; D1 = A; D2 = B;
elseif m <= stop
    [T,D] = eig(A,B);
    for i = 1:m
        T(:,i) = T(:,i)/norm(T(:,i),2);
    end
    D1 = D;
    D2 = eye(m);
else
    % ---------------
    % Set parameters
    % ---------------
    zeta = 2*floor(log2(s)+1);
    eta = min([4*pi*beta*epsilon.^2/(315*sqrt(8)*omega*n.^alpha); 1/(2*(log(n)/log(5/4)))]);
    delta = min([sqrt(theta/10)*epsilon.^2/(7200*n.^(2*alpha+3)); theta/(2*(theta + 10*n.^6*zeta)); sqrt(theta/10)*eta.^2/(288*n.^(2*alpha+3))]);
    p = ceil(log2(n/epsilon)); % p = ceil(max([7; log2((105*n.^alpha/epsilon)-1); -2*log2(-(1/2)*log2(1-(epsilon/(105*n.^alpha)))); 1+log2(log2(delta*pi*epsilon/(12*n.^alpha*m*omega+delta*pi*epsilon))/log2(1-(epsilon/(105*n.^alpha))))]));
    % ---------------------------------------------------------------
    % Start searching over vertical grid lines (via a binary search)
    % ---------------------------------------------------------------
    index = max(ceil(s_v/2),1);
    left_index = 1;
    right_index = s_v;
    lines_checked = 0;
    split_failed = false;
    line_found = false;
    no_vertical_split = false;
    if index == 1
        no_vertical_split = true;
    end
    while line_found == false && no_vertical_split == false
        h = g(1)+(index-1)*omega;
        Ascript = A-(h-1)*B;
        Bscript = A-(h+1)*B;
        [A_p,B_p] = irs(Ascript,Bscript,p);
        [~,R1,R2,~] = grurv(A_p+B_p,A_p,-1,1);
        values = zeros(m,1);
        for i=1:m
            values(i) = abs(R2(i,i)/R1(i,i));
        end
        lines_checked = lines_checked + 1;
        k = nnz(values(:) >= sqrt(theta/(10*zeta))*(1-delta)); %sqrt(theta/(10*zeta))*(1-delta)/n.^3);
        if k < m/5 % search to the left
            right_index = index;
            index = max(floor((left_index + index)/2),1);
        elseif k > 4*m/5 % search to the right
            left_index = index;
            index = min(floor((right_index + index)/2),s_v);   
        else
           line_found = true;
        end
        if right_index-left_index <= 1 || lines_checked >= floor(log2(s_v)+1)
            no_vertical_split = true;
            vertical_lines_checked = lines_checked;
        end
    end
    % ----------------------------------
    % Apply vertical split if it exists
    % ----------------------------------
    if line_found == true
        disp('vertical split')
        if m > 3 
            splits = [splits; k/m];
        end
        [U_R1,U_L1] = deflate(Ascript,Bscript,p,k);
        Ascript = A-(h+1)*B;
        Bscript = A-(h-1)*B;
        [U_R2,U_L2] = deflate(Ascript,Bscript,p,m-k);
        A_11 = U_L1'*A*U_R1;
        B_11 = U_L1'*B*U_R1;
        A_22 = U_L2'*A*U_R2;
        B_22 = U_L2'*B*U_R2;
        g_R = [h;g(2);g(3);g(4)];
        g_L = [g(1);h;g(3);g(4)];
        flop_count = flop_count + lines_checked*m.^3;
        [T_R,D_R1,D_R2,splits,fails,flop_count] = dnc_eig(n,A_11,B_11,4*epsilon/5,alpha,g_R,beta/3,omega,theta,splits,fails,flop_count,stop);
        [T_L,D_L1,D_L2,splits,fails,flop_count] = dnc_eig(n,A_22,B_22,4*epsilon/5,alpha,g_L,beta/3,omega,theta,splits,fails,flop_count,stop);
        T = [U_R1*T_R U_R2*T_L];
        for i = 1:m
            T(:,i) = T(:,i)/norm(T(:,i),2);
        end
        D1 = [D_R1 zeros(k,m-k); zeros(m-k,k) D_L1];
        D2 = [D_R2 zeros(k,m-k); zeros(m-k,k) D_L2];
    else
        % ------------------------------------------------------------
        % Switch to searching over horizontal grid lines if necessary
        % ------------------------------------------------------------
        lines_checked = 0;
        index = max(ceil(s_h/2),1);
        bottom_index = 1;
        top_index = s_h;
        while line_found == false && split_failed == false
            h = g(3)+(index-1)*omega;
            Ascript = A-1i*(h-1)*B;
            Bscript = A-1i*(h+1)*B;
            [A_p,B_p] = irs(Ascript,Bscript,p);
            [~,R1,R2,~] = grurv(A_p+B_p,A_p,-1,1);
            for i = 1:m
                values(i) = abs(R2(i,i)/R1(i,i));
            end
            lines_checked = lines_checked+1;
            k = nnz(values(:) >= sqrt(theta/(10*zeta))*(1-delta)); %sqrt(theta/(10*zeta))*(1-delta)/n.^3);
            if k<m/5 % search below
                top_index = index;
                index = max(floor((index+bottom_index)/2),1);
            elseif k>4*m/5 % search above
                bottom_index = index;
                index = min(floor((index+top_index)/2),s_h);
            else
                line_found = true;
            end
            if top_index-bottom_index <= 1 || lines_checked >= floor(log2(s_h)+1)
                disp('error: no dividing line found')
                fails = [fails; m];
                split_failed = true;
            end
        end
        % ---------------------
        % Call QZ if necessary
        % ---------------------
        if split_failed 
            [T,D1] = eig(A,B);
            D2 = eye(m);
            flop_count = flop_count + (lines_checked+vertical_lines_checked)*m.^3;
        else
            % -----------------------
            % Apply horizontal split
            % -----------------------
            disp('horizontal split')
            if m > 3
                splits = [splits; k/m];
            end
            [U_R1,U_L1] = deflate(Ascript,Bscript,p,k);
            Ascript = A-1i*(h+1)*B;
            Bscript = A-1i*(h-1)*B;
            [U_R2,U_L2] = deflate(Ascript,Bscript,p,m-k);
            A_11 = U_L1'*A*U_R1;
            B_11 = U_L1'*B*U_R1;
            A_22 = U_L2'*A*U_R2;
            B_22 = U_L2'*B*U_R2;
            g_T = [g(1);g(2);h;g(4)];
            g_B = [g(1);g(2);g(3);h];
            flop_count = flop_count + (lines_checked+vertical_lines_checked)*m.^3;
            [T_A,DA1,DA2,splits,fails,flop_count] = dnc_eig(n,A_11,B_11,4*epsilon/5,alpha,g_T,beta/3,omega,theta,splits,fails,flop_count,stop);
            [T_B,DB1,DB2,splits,fails,flop_count] = dnc_eig(n,A_22,B_22,4*epsilon/5,alpha,g_B,beta/3,omega,theta,splits,fails,flop_count,stop);
            T = [U_R1*T_A U_R2*T_B];
            for i = 1:m
                T(:,i) = T(:,i)/norm(T(:,i),2);
            end
            D1 = [DA1 zeros(k,m-k); zeros(m-k,k) DB1];
            D2 = [DA2 zeros(k,m-k); zeros(m-k,k) DB2];
        end
    end
end

