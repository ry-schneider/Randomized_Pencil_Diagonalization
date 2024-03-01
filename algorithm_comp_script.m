% ------------------------------------------------------------------------
% This script applies ALGORITHM_COMP to a series of randomly drawn pencils
% (A,B). In each case, B is modified to be singular before applying either
% divide-and-conquer algorithm. 
%
% Two values for the backward diagaonlization error are used: 10^{-5} and 
% 10^{-10}. Histograms provide performance data for both algorithms.
% -------------------------------------------------------------------------
n = 1000; % problem size
samples = 10; % number of times each pencil is run through ALGORITHM_COMP
trials = 20; % number of random draws of (A,B)
%-------------------------
% Initialize error arrays
%-------------------------
eig_errors = [];
diag_errors = [];
flop_array = [];
banks_fails = [];
banks_splits = [];
inverse_free_fails = [];
inverse_free_splits = [];
condition_numbers = [];
qz_errors = [];
qr_errors = [];
eig_errors2 = [];
diag_errors2 = [];
flop_array2 = [];
banks_fails2 = [];
banks_splits2 = [];
inverse_free_fails2 = [];
inverse_free_splits2 = [];
condition_numbers2 = [];
qz_errors2 = [];
qr_errors2 = [];
%--------------
% Start trials
%--------------
for i = 1:trials
    disp(i);
    %--------------------
    % Draw random inputs
    %--------------------
    A = randn(n)+1i*randn(n);
    B = randn(n)+1i*randn(n);
    %-------------------------
    % Modify B to be singular
    %-------------------------
    [U,Sigma,V] = svd(B);
    B = B-U(:,n)*Sigma(n,n)*V(:,n)';
    %---------------------------------------------------------------
    % Run both versions of divide-and-conquer with error = 0.000001
    %---------------------------------------------------------------
    for j = 1:samples
        [e_1,d_1,fails1,splits1,flops1,cond_1,qz_error,e_2,d_2,fails2,splits2,flops2,cond_2,qr_error] = algorithm_comp(A,B,0.00001,250);
        if size(fails1,1) == 0 && size(fails2,1) == 0
            eig_errors = [eig_errors; log10(e_1) log10(e_2)];
            diag_errors = [diag_errors; log10(d_1) log10(d_2)];
            inverse_free_splits = [inverse_free_splits; splits1];
            banks_splits = [banks_splits; splits2];
            flop_array = [flop_array; flops1/1250000000 flops2/1250000000];
            condition_numbers = [condition_numbers; log10(cond_1) log10(cond_2)];
            qz_errors = [qz_errors; qz_error];
            qr_errors = [qr_errors; qr_error];
        else
            inverse_free_fails = [inverse_free_fails; fails1];
            banks_fails = [banks_fails; fails2];
        end
    end
    %-------------------------------------------------------------------
    % Run both versions of divide-and-conquer with error = 0.0000000001
    %-------------------------------------------------------------------
    for j = 1:samples
        [e_1,d_1,fails1,splits1,flops1,cond_1,qz_error,e_2,d_2,fails2,splits2,flops2,cond_2,qr_error] = algorithm_comp(A,B,0.0000000001,250);
        if size(fails1,1) == 0 && size(fails2,1) == 0
            eig_errors2 = [eig_errors2; log10(e_1) log10(e_2)];
            diag_errors2 = [diag_errors2; log10(d_1) log10(d_2)];
            inverse_free_splits2 = [inverse_free_splits2; splits1];
            banks_splits2 = [banks_splits2; splits2];
            flop_array2 = [flop_array2; flops1/1250000000 flops2/1250000000];
            condition_numbers2 = [condition_numbers2; log10(cond_1) log10(cond_2)];
            qz_errors2 = [qz_errors2; qz_error];
            qr_errors2 = [qr_errors2; qr_error];
        else
            inverse_free_fails2 = [inverse_free_fails2; fails1];
            banks_fails2 = [banks_fails2; fails2];
        end
    end
end
%------------------------
% Plot error histograms
%------------------------
figure
tiledlayout(2,4,'TileSpacing','loose');
nexttile
histogram(diag_errors(:,1),'FaceColor',[0 0.4470 0.7410])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
% ------------------------------------------------------------------------
nexttile
histogram(diag_errors(:,2),'FaceColor',[0 0.4470 0.7410])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
% ------------------------------------------------------------------------
nexttile
histogram(diag_errors2(:,1),'FaceColor',[0 0.4470 0.7410])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
% ------------------------------------------------------------------------
nexttile
histogram(diag_errors2(:,2),'FaceColor',[0 0.4470 0.7410])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
% ------------------------------------------------------------------------
nexttile
val = mean(log10(qz_errors));
histogram(eig_errors(:,1),'FaceColor',[0.4660 0.6740 0.1880])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
xline(val,'--r','QZ','Interpreter','latex','LineWidth',2,'FontSize',16)
% ------------------------------------------------------------------------
nexttile
val = mean(log10(qr_errors));
histogram(eig_errors(:,2),'FaceColor',[0.4660 0.6740 0.1880])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
xline(val,'--r','QR','Interpreter','latex','LineWidth',2,'FontSize',16)
% ------------------------------------------------------------------------
nexttile
val = mean(log10(qz_errors2));
histogram(eig_errors2(:,1),'FaceColor',[0.4660 0.6740 0.1880])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
xline(val,'--r','QZ','Interpreter','latex','LineWidth',2,'FontSize',16)
% ------------------------------------------------------------------------
nexttile
val = mean(log10(qr_errors2));
histogram(eig_errors2(:,2),'FaceColor',[0.4660 0.6740 0.1880])
set(gca,'XTick',[-0.6 -0.2 0.2],'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','Interpreter','latex','FontSize',18)
xline(val,'--r','QR','Interpreter','latex','LineWidth',2,'FontSize',16)
%--------------------------------------------------------------------------
figure 
tiledlayout(3,4)
nexttile
histogram(inverse_free_splits,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$k/m$','Interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
split_size = size(inverse_free_splits,1);
formatSpec = "%d";
str = sprintf(formatSpec,split_size);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(banks_splits,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$k/m$','Interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
split_size = size(banks_splits,1);
formatSpec = "%d";
str = sprintf(formatSpec,split_size);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(inverse_free_splits2,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$k/m$','Interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
split_size = size(inverse_free_splits2,1);
formatSpec = "%d";
str = sprintf(formatSpec,split_size);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(banks_splits2,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$k/m$','Interpreter','latex','FontSize',14)
xL=xlim;
yL=ylim;
split_size = size(banks_splits2,1);
formatSpec = "%d";
str = sprintf(formatSpec,split_size);
text(0.98*xL(2),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(flop_array(:,1),'FaceColor',[0.6350 0.0780 0.1840])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('Lines Checked','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(flop_array(:,2),'FaceColor',[0.6350 0.0780 0.1840])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('Lines Checked','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(flop_array2(:,1),'FaceColor',[0.6350 0.0780 0.1840])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('Lines Checked','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(flop_array2(:,2),'FaceColor',[0.6350 0.0780 0.1840])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('Lines Checked','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(condition_numbers(:,1),'FaceColor',[0.4940 0.1840 0.5560])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$\log_{10}(\kappa_2(T))$','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(condition_numbers(:,2),'FaceColor',[0.4940 0.1840 0.5560])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$\log_{10}(\kappa_2(T))$','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(condition_numbers2(:,1),'FaceColor',[0.4940 0.1840 0.5560])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$\log_{10}(\kappa_2(T))$','Interpreter','latex','FontSize',14)
% ------------------------------------------------------------------------
nexttile
histogram(condition_numbers2(:,2),'FaceColor',[0.4940 0.1840 0.5560])
set(gca,'TickLabelInterpreter','latex','FontSize',14)
xlabel('$\log_{10}(\kappa_2(T))$','Interpreter','latex','FontSize',14)
