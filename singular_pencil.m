% -----------------------------------------------------------------------
% This script applies RPD to a 4 x 4 singular pencil (A,B) taken from a
% paper of Lotz and Noferini (Foundations of Computational Math 2020). Note
% that (A,B) has only one eigenvalue at lambda = 1.
%
% Each run of RPD uses error = 10^-6, and we again produce histograms of
% performance data for the algorithm. 
%
% samples is the number of runs of RPD (default is 500).
% -----------------------------------------------------------------------
format long
% ----------------------------------------
% Construct singular pencil and normalize
% ----------------------------------------
A = [2 -1 -5 -1; 6 -2 -11 -2; 5 0 -2 0; 3 1 3 1];
B = -[-1 1 4 2; -2 3 12 6; 1 3 11 6; 2 2 7 4];
c = max(norm(A,2), norm(B,2));
Anorm = A/c; Bnorm = B/c;
samples = 500;
split_record = [];
fail_record = [];
diag_error = zeros(samples,1);
eig_error = zeros(samples,1);
errors = zeros(4,1);
% --------------------------------------------------
% For each sample, call RPD and record error/splits
% --------------------------------------------------
for i = 1:samples
    [S,T,D,splits,fails] = rpd(Anorm,Bnorm,0.000001);
    split_record = [split_record; splits];
    fail_record = [fail_record; fails];
    diag_error(i) = log10(max(norm(Anorm-S*D/T,2), norm(Bnorm-S/T,2)));
    eigs = diag(D);
    for j = 1:4
        errors(j) = log10(abs(eigs(j)-1));
    end
    eig_error(i) = min(errors);
end
% ------------------------------------------------------------------
% Create histograms of splits and eigenvalue/diagonalization errors
% ------------------------------------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(eig_error,'FaceColor',[0.4660 0.6740 0.1880])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}\left( \min_j|\lambda_j - 1| \right)$','interpreter','latex','FontSize',18)
title('Eigenvalue Accuracy','interpreter','latex','FontSize',18)
% ------------------------------------------------------------------------------------
nexttile
histogram(diag_error)
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$\log_{10}$(error)','interpreter','latex','FontSize',18)
title('Diagonalization Accuracy','interpreter','latex','FontSize',18)
% xline(-6,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LabelHorizontalAlignment','left','LineWidth',2,'FontSize',14)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(diag_error(:) > -6));
text(xL(2)-0.02*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% ------------------------------------------------------------------------------------
nexttile
histogram(split_record,'FaceColor',[0.8500 0.3250 0.0980])
set(gca,'TickLabelInterpreter','latex','FontSize',18)
xlabel('$k/m$','interpreter','latex','FontSize',18)
title('Split Size','interpreter','latex','FontSize',18)



