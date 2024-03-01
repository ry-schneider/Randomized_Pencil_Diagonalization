% -----------------------------------------------------------------------
% This script applies RPD to a pencil (J,I) where J is a Jordan block with
% one eigenvalue at zero. 
%
% RPD is run with three different values of error, and plots of 
% performance data are made.
%
% n is the size of the problem (default is 50 x 50).
% samples is the number of runs of RPD (default is 500).
% -----------------------------------------------------------------------
n = 50;
samples = 500;
% -----------------------------
% Construct the Jordan block J
% -----------------------------
J = zeros(n);
for i = 1:n-1
    J(i,i+1) = 1;
end
diag_error = zeros(samples,3);
flop_count = zeros(samples,3);
eigs = zeros(n,3);
splits1 = [];
splits2 = [];
splits3 = [];
fails1 = [];
fails2 = [];
fails3 = [];
% ------------------------------------------------------------------------
% Repeatedly apply RPD, recording diagonalization error, split sizes, and
% calls to QZ
% ------------------------------------------------------------------------
for i = 1:samples
    [S1,T1,D1,splits,fails,flops] = rpd(J,eye(n),0.01);
    splits1 = [splits1; splits];
    fails1 = [fails1; fails];
    flop_count(i,1) = flops/166732;
    [S2,T2,D2,splits,fails,flops] = rpd(J,eye(n),0.001);
    splits2 = [splits2; splits];
    fails2 = [fails2; fails];
    flop_count(i,2) = flops/166732;
    [S3,T3,D3,splits,fails,flops] = rpd(J,eye(n),0.0001);
    splits3 = [splits3; splits];
    fails3 = [fails3; fails];
    flop_count(i,3) = flops/166732;
    diag_error(i,1) = log10(max(norm(J-S1*D1/T1,2),norm(eye(n)-S1/T1,2)));
    % ---------------------------------------------------------------
    % If the diagonalization is successful, record the corresponding
    % eigenvalue approximations.
    % ---------------------------------------------------------------
    if diag_error(i,1) < -2 
        eigs(:,1) = diag(D1);
    end
    diag_error(i,2) = log10(max(norm(J-S2*D2/T2,2),norm(eye(n)-S2/T2,2)));
    if diag_error(i,2) < -3 
        eigs(:,2) = diag(D2);
    end
    diag_error(i,3) = log10(max(norm(J-S3*D3/T3,2),norm(eye(n)-S3/T3,2)));
    if diag_error(i,3) < -4
        eigs(:,3) = diag(D3);
    end
end
% -----------------------------
% Plot approximate eigenvalues
% -----------------------------
figure
tiledlayout(1,3)
nexttile
scatter(real(eigs(:,1)),imag(eigs(:,1)),'filled');
axis([-1 1 -1 1]);
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
xlabel('Real','interpreter','latex','FontSize',18);
ylabel('Imaginary','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
% --------------------------------------------------------------------------------
nexttile
scatter(real(eigs(:,2)),imag(eigs(:,2)),'filled');
axis([-1 1 -1 1]);
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
xlabel('Real','interpreter','latex','FontSize',18);
ylabel('Imaginary','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
% --------------------------------------------------------------------------------
nexttile
scatter(real(eigs(:,3)),imag(eigs(:,3)),'filled');
axis([-1 1 -1 1]);
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
xlabel('Real','interpreter','latex','FontSize',18);
ylabel('Imaginary','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
% -------------------------------------------
% Create histograms of diagonalization error
% -------------------------------------------
figure
tiledlayout(1,3)
nexttile
histogram(diag_error(:,1))
title('$\varepsilon = 10^{-2}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
xline(-2,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(diag_error(:,1) > -2));
text(xL(2)-0.025*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% --------------------------------------------------------------------------------
nexttile
histogram(diag_error(:,2))
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
xline(-3,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(diag_error(:,2) > -3));
text(xL(2)-0.025*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
% --------------------------------------------------------------------------------
nexttile
histogram(diag_error(:,3))
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('$\log_{10}$(diag error)','interpreter','latex','FontSize',18);
% xline(-4,'--r','$\log_{10}(\varepsilon)$','Interpreter','latex','LineWidth',2,'LabelHorizontalAlignment','left','FontSize',16)
xL=xlim;
yL=ylim;
formatSpec = "Fails: %d";
str = sprintf(formatSpec,nnz(diag_error(:,3) > -4));
text(xL(2)-0.025*(xL(2)-xL(1)),0.98*yL(2),str,'interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top','FontSize',16)
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
% --------------------------------------------------------------------------------
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
% --------------------------------------------------------------------------------
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
% --------------------------------------------------------------------------------
nexttile
histogram(flop_count(:,2),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\varepsilon = 10^{-3}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',18);
% ylabel('Frequency','interpreter','latex','FontSize',14);
% --------------------------------------------------------------------------------
nexttile
histogram(flop_count(:,3),'FaceColor',[0.6350 0.0780 0.1840]);
title('$\varepsilon = 10^{-4}$','interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','FontSize',18);
xlabel('Relative Efficiency Factor','interpreter','latex','FontSize',18);
