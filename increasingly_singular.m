% -----------------------------------------------------------------------
% This script applies RPD to a series of pencils (A,B_tau) where B_tau
% becomes singular as tau approaches one. Initially (i.e., for tau = 0) 
% the pencil is the planted spectrum example.
%
% We first apply RPD to produce eigenvalue approxiamtions for twenty one
% values of tau. We then apply RPD repeatedly to compute empirical
% statistics for the algorithm's performance with three values of tau (0, 
% 1/2, and 1). In all runs, error = 10^-4.
%
% n is the size of the problem (default is 100 x 100).
% samples is the number of runs of RPD used to construct performance
% histograms (default is 1000).
% -----------------------------------------------------------------------
n = 100;
% -----------------------------------------------------
% Construct (A,B) to have evenly spaced eigs on [-2,2]
% -----------------------------------------------------
d = zeros(n,1);
for i = 1:n
    d(i) = -2+(i-1)*4/(n-1);
end
E = diag(d);
X = 1/sqrt(2)*(randn(n)+1i*randn(n));
Y = 1/sqrt(2)*(randn(n)+1i*randn(n));
A = X*E/Y;
B = X/Y;
% -----------------
% Compute SVD of B
% -----------------
[U,Sigma,V] = svd(B);
Z = zeros(n);
Z(n,n) = Sigma(n,n);
box_plot_errors = zeros(n,21);
box_plot_grouping = zeros(21,1);
diag_error = ones(21,1);
% ---------------------------------------------------------------------
% Find eigenvalue approximations for tau = (j-1)/20 and compute errors
% ---------------------------------------------------------------------
for j = 1:21
    Atest = A;
    Btest = B - ((j-1)/20)*U*Z*V';
    c = max(norm(Atest,2),norm(Btest,2));
    Atest = Atest/c; 
    Btest = Btest/c;
    while diag_error(j) > 0.0001
        [S,T,D,~,~] = rpd(Atest,Btest,0.0001);
        diag_error(j) = max(norm(Atest-S*D/T,2), norm(Btest-S/T,2));
    end
    true_eigs = eig(Atest,Btest);
    eig_approx = diag(D);
    true_eigs = sort(true_eigs,'ComparisonMethod','real');
    eig_approx = sort(eig_approx,'ComparisonMethod','real');
    for i = 1:n
        box_plot_errors(i,j) = log10(abs(true_eigs(i)-eig_approx(i)));
    end 
    box_plot_grouping(j) = (j-1)/20;
end
samples = 1000;
tau_num = 3;
split_sizes = zeros(tau_num,1);
fail_sizes = zeros(tau_num,1);
d_errors = zeros(samples,tau_num);
flop_count = zeros(samples,tau_num);
split_record = [];
fail_record = [];
% ------------------------------------------------------------------------
% Now repeatedly apply RPD for three different values of tau, each time
% recording diagonalization error, split sizes, and calls to QZ
% ------------------------------------------------------------------------
for j = 1:tau_num
    Atest = A;
    Btest = B - (j-1)/(tau_num-1)*U*Z*V';
    c = max(norm(Atest,2),norm(Btest,2));
    Atest = Atest/c; Btest = Btest/c;
    for k = 1:samples
        [S,T,D,splits,fails,flops] = rpd(Atest,Btest,0.0001);
        split_record = [split_record; splits];
        fail_record = [fail_record; fails];
        d_errors(k,j) = log10(max(norm(Atest-S*D/T,2),norm(Btest-S/T,2)));
        flop_count(k,j) = flops/1333464;
    end
    split_sizes(j) = size(split_record,1);
    fail_sizes(j) = size(fail_record,1);
end
% --------------------------------------
% Create box plots of eigenvalue errors
% --------------------------------------
figure
nexttile
boxplot(box_plot_errors,box_plot_grouping)
xlabel('$\tau$','interpreter','latex','FontSize',14)
ylabel('$\log_{10}$(error)','interpreter','latex','FontSize',14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)
% -------------------------------------------
% Create histograms of diagonalization error
% -------------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(d_errors(:,1))
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 0$','interpreter','latex','FontSize',14)
xlabel('$\log_{10}\left( \max \left\{ ||A-SDT^{-1}||_2, ||B_{\tau}-SIT^{-1}||_2 \right\} \right)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xline(-4,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',14)
nexttile
histogram(d_errors(:,2))
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 0.5$','interpreter','latex','FontSize',14)
xlabel('$\log_{10}\left( \max \left\{ ||A-SDT^{-1}||_2, ||B_{\tau}-SIT^{-1}||_2 \right\} \right)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xline(-4,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',14)
nexttile
histogram(d_errors(:,3))
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 1$','interpreter','latex','FontSize',14)
xlabel('$\log_{10}\left( \max \left\{ ||A-SDT^{-1}||_2, ||B_{\tau}-SIT^{-1}||_2 \right\} \right)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xline(-4,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',14)
% ---------------------------------
% Create histograms of split sizes
% ---------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(split_record(1:split_sizes(1)),'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 0$','interpreter','latex','FontSize',14)
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
str = sprintf(formatSpec,split_sizes(1));
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
nexttile
histogram(split_record(split_sizes(1)+1:split_sizes(2)),'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 0.5$','interpreter','latex','FontSize',14)
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
split_num = split_sizes(2)-split_sizes(1);
str = sprintf(formatSpec,split_num);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
nexttile
histogram(split_record(split_sizes(2)+1:split_sizes(3)),'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('$\tau = 1$','interpreter','latex','FontSize',14)
xlabel('Split Size $(k/m)$','interpreter','latex','FontSize',14)
ylabel('Frequency','interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
formatSpec = "Total: %d";
split_num = split_sizes(3)-split_sizes(2);
str = sprintf(formatSpec,split_num);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
% ----------------------------------------
% Create histograms of pseudo-flop counts
% ----------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(flop_count(:,1),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\tau = 0$','interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',14);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14);
text()
nexttile
histogram(flop_count(:,2),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\tau = 0.5$','interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',14);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14);
nexttile
histogram(flop_count(:,3),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\tau = 1$','interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize',14);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14);
