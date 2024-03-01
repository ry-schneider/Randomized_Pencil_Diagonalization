% ------------------------------------------------------------------------
% This script applies RPD to a planted spectrum example. Here, (A,B) is
% constructed to have evenly spaced (real) eigenvalues in [-2,2].
%
% RPD is run with three different values of error, and plots of 
% performance data are made.
%
% n is the size of the problem (default is 50 x 50).
% samples is the number of runs of RPD (default is 500)
% ------------------------------------------------------------------------
n = 50;
samples = 500;
errors = zeros(samples,3);
eig_errors = zeros(n,3);
flop_count = zeros(samples,3);
% ----------------
% Construct (A,B)
% ----------------
d = zeros(n,1);
for i = 1:n
    d(i) = -2+(i-1)*4/(n-1);
end
D = diag(d);
X = 1/sqrt(2)*(randn(n)+1i*randn(n));
Y = 1/sqrt(2)*(randn(n)+1i*randn(n));
A = X*D/Y;
B = X/Y;
% -----------------------------
% Normalize before call to rpd
% -----------------------------
c = max(norm(A,2),norm(B,2));
A = A/c; B = B/c;
splits1 = [];
splits2 = [];
splits3 = [];
fails1 = [];
fails2 = [];
fails3 = [];
% ---------------------------------------------------------------------
% For each sample, find diagonalization, compute residuals, and record
% split sizes/calls to QZ
% ---------------------------------------------------------------------
for i = 1:samples
    [S1,T1,D1,splits,fails,flops] = rpd(A,B,0.01);
    splits1 = [splits1; splits];
    fails1 = [fails1;fails];
    flop_count(i,1) = flops/166732;
    [S2,T2,D2,splits,fails,flops] = rpd(A,B,0.001);
    splits2 = [splits2; splits];
    fails2 = [fails2; fails];
    flop_count(i,2) = flops/166732;
    [S3,T3,D3,splits,fails,flops] = rpd(A,B,0.0001);
    splits3 = [splits3; splits];
    fails3 = [fails3; fails];
    flop_count(i,3) = flops/166732;
    errors(i,1) = max(norm(A-S1*D1/T1,2), norm(B-S1*eye(n)/T1,2));
    % ---------------------------------------------------------------
    % If the diagonalization is successful, record the corresponding
    % eigenvalue approximations.
    % ---------------------------------------------------------------
    if errors(i,1) < 0.01
        d1 = diag(D1);
    end
    errors(i,2) = max(norm(A-S2*D2/T2,2), norm(B-S2*eye(n)/T2,2));
    if errors(i,2) < 0.001 
        d2 = diag(D2);
    end
    errors(i,3) = max(norm(A-S3*D3/T3,2), norm(B-S3*eye(n)/T3,2));
    if errors(i,3) < 0.0001
        d3 = diag(D3);
    end
end
% ----------------------------------------
% Compute eigenvalue approximation errors
% ----------------------------------------
d1 = sort(d1,'ComparisonMethod','real');
d2 = sort(d2,'ComparisonMethod','real');
d3 = sort(d3,'ComparisonMethod','real');
indices = zeros(n,1);
for j = 1:n
    eig_errors(j,1) = d1(j) - d(j);
    eig_errors(j,2) = d2(j) - d(j);
    eig_errors(j,3) = d3(j) - d(j);
    indices(j) = j;
end
% --------------------------
% Set colors for each error
% --------------------------
c1 = zeros(n,1);
c2 = zeros(n,1);
c3 = zeros(n,1);
for k = 1:n
    c1(k) = log10(abs(eig_errors(k,1)));
    c2(k) = log10(abs(eig_errors(k,2)));
    c3(k) = log10(abs(eig_errors(k,3)));
end
% -----------------------------
% Plot approximate eigenvalues
% -----------------------------
figure
tiledlayout(1,3)
nexttile
scatter(real(d1),imag(d1),[],c1,'filled')
colorbar
xlim([-2 2]);
ylim([-0.3 0.3]);
clim([-5,1]);
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
xlabel('Real', 'interpreter', 'latex','FontSize',18)
ylabel('Imaginary', 'interpreter','latex','FontSize',18)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xticks([-2 0 2])
yticks([-0.3 0 0.3])
set(colorbar,'YTick',[-5 -3 -1 1],'TickLabelInterpreter','latex','FontSize',18);
% -------------------------------------------------------------------------------------
nexttile
scatter(real(d2),imag(d2),[],c2,'filled')
colorbar
xlim([-2 2]);
ylim([-0.3 0.3]);
clim([-5 1]);
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
xlabel('Real', 'interpreter', 'latex','FontSize',18)
ylabel('Imaginary', 'interpreter','latex','FontSize',18)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xticks([-2 0 2])
yticks([-0.3 0 0.3])
set(colorbar,'YTick',[-5 -3 -1 1],'TickLabelInterpreter','latex','FontSize',18);
% -------------------------------------------------------------------------------------
nexttile
scatter(real(d3),imag(d3),[],c3,'filled')
colorbar
xlim([-2 2]);
ylim([-0.3 0.3]);
clim([-5 1]);
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
xlabel('Real', 'interpreter', 'latex','FontSize',18)
ylabel('Imaginary', 'interpreter','latex','FontSize',18)
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xticks([-2 0 2])
yticks([-0.3 0 0.3])
set(colorbar,'YTick',[-5 -3 -1 1],'TickLabelInterpreter','latex','FontSize',18);
% -------------------------------------------
% Create histograms of diagonalization error
% -------------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(log10(errors(:,1)))
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
xline(-2,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(log10(errors(:,1)) > -2));
text(xL(2)-0.025*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% -------------------------------------------------------------------------------------
nexttile
histogram(log10(errors(:,2)))
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
xline(-3,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(log10(errors(:,2)) > -3));
text(xL(2)-0.025*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% -------------------------------------------------------------------------------------
nexttile
histogram(log10(errors(:,3)))
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
xline(-4,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(log10(errors(:,3)) > -4));
text(xL(2)-0.02*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% ---------------------------------
% Create histograms of split sizes
% ---------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(splits1,'FaceColor',[0.8500 0.3250 0.0980])
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',18);
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
str = sprintf(formatSpec,size(splits1,1));
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% -------------------------------------------------------------------------------------
nexttile
histogram(splits2,'FaceColor',[0.8500 0.3250 0.0980])
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',18);
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
str = sprintf(formatSpec,size(splits2,1));
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% -------------------------------------------------------------------------------------
nexttile
histogram(splits3,'FaceColor',[0.8500 0.3250 0.0980]);
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',18);
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
str = sprintf(formatSpec,size(splits3,1));
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% ----------------------------------------
% Create histograms of pseudo-flop counts
% ----------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(flop_count(:,1),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',18);
% -------------------------------------------------------------------------------------
nexttile
histogram(flop_count(:,2),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',18);
% -------------------------------------------------------------------------------------
nexttile
histogram(flop_count(:,3),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',18);
