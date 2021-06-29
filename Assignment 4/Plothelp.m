%% Plotting Results
close all
%% Single Plot
close all
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.1, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])


plot((avgerr1), 'Color', '#D95F69', 'Linewidth', 2);
plot((err1),'-.','Color', '#6ABABD', 'Linewidth', 0.1);
plot((err2),'-.','Color', '#F9C139', 'Linewidth', 0.1);
plot((err3),'-.','Color', '#75D606', 'Linewidth', 0.1);
plot((err4),'-.','Color', '#E67F5A', 'Linewidth', 0.1);
% yline((1e-3), 'r--')


plot((avgerr), 'Color', '#D95F69', 'Linewidth', 2);
% plot((avgerravg),'k--', 'Linewidth', 0.1);

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||\bar x_f^i - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

% leg = legend('$$\mathbf{\varphi = 0}$$','$$\mathbf{\varphi = 1}$$','$$\mathbf{\varphi = 5}$$','$$\mathbf{\varphi = 100}$$','$$\mathbf{\varphi = \infty}$$', 'Interpreter', 'Latex');
leg = legend('$$\mathbf{\alpha = 0.1}$$','$$\mathbf{\alpha = 0.05}$$','$$\mathbf{\alpha = 0.01}$$','$$\mathbf{\alpha = 0.005}$$', 'Interpreter', 'Latex');
leg.FontSize=20;

tit = title('\textbf{Average Error Evolution Normalized Against Centralized Solution, Consenus Iterations, $$\varphi = 20$$.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xf1v)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.

%% Single Plot
close all
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])

plot((vecnorm(xf1v-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2); 
% plot((vecnorm(xf3v-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2); 
% plot((vecnorm(xf4v-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2); 
% % plot((vecnorm(xfv-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||x_f - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

% leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');
leg = legend('$$\rho = 1$$', '$$\rho = 4$$', '$$\rho = 7$$', '$$\rho = 12$$','Interpreter', 'Latex');

leg.FontSize=30;

tit = title('\textbf{Error Evolution Normalized Against Centralized Solution.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xf1v)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.
%% Two Plots Side by Side

figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.025, 0.1, 0.95, 0.9]);
sp1 = subplot(1,2,1);
plot((vecnorm(xf1v_05-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v_05-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2); 
plot((vecnorm(xf3v_05-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2); 
plot((vecnorm(xf4v_05-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2); 
xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel(['$$\mathbf{\frac{||x_f^i - x_f^*||_2}{||x_f^*||}}$$'], 'Interpreter', 'Latex', 'Fontsize', 60);

leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');
leg.FontSize=30;

xlim([0,length(xf1v_05)])
grid on; grid minor;
set(gca, 'Fontsize', 16, 'Yscale', 'log')
title(['$$\alpha = \frac{k-1}{k+2}$$.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')

sp2 = subplot(1,2,2);
plot((vecnorm(xf1v-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2); 
plot((vecnorm(xf3v-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2); 
plot((vecnorm(xf4v-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2); 
xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');
leg.FontSize=30;

xlim([0,length(xf1v)])
grid on; grid minor;
set(gca, 'Fontsize', 16, 'Yscale', 'log')
title(['$$\alpha = 0.75\cdot\frac{k-1}{k+2}$$.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')

sgtitle('\textbf{Error Evolution Normalized w.r.t. Centralized Solution. Gradient not Normalized, Four Constraints.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');

sp1.Position = [0.0814 0.1057 0.4239 0.7432];
sp2.Position = [0.5548 0.1057 0.4239 0.7432];


%% Four Plots In a Grid

figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.025, 0.1, 0.95, 0.85]);
sp1 = subplot(2,2,1);
plot(0:5,x1(1:4:end-3), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot(0:5,x2(1:4:end-3), 'Color', '#6593A6', 'Linewidth', 2); 
plot(0:5,x3(1:4:end-3), 'Color', '#A9BF04', 'Linewidth', 2); 
plot(0:5,x4(1:4:end-3), 'Color', '#D95F69', 'Linewidth', 2); 

ylabel(['$$\mathbf{x^i(t)}$$'], 'Interpreter', 'Latex', 'Fontsize', 60);

leg = legend('$$\mathbf{x^1(t)}$$','$$\mathbf{x^2(t)}$$','$$\mathbf{x^3(t)}$$','$$\mathbf{x^4(t)}$$', 'Interpreter', 'Latex');
leg.FontSize=15; leg.Orientation='Horizontal';


xlim([0, 5])
grid on; grid minor;
set(gca, 'Fontsize', 16)
title(['Component 1.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')

sp2 = subplot(2,2,2);
plot(0:5,x1(2:4:end-2), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot(0:5,x2(2:4:end-2), 'Color', '#6593A6', 'Linewidth', 2); 
plot(0:5,x3(2:4:end-2), 'Color', '#A9BF04', 'Linewidth', 2); 
plot(0:5,x4(2:4:end-2), 'Color', '#D95F69', 'Linewidth', 2); 

leg = legend('$$\mathbf{x^1(t)}$$','$$\mathbf{x^2(t)}$$','$$\mathbf{x^3(t)}$$','$$\mathbf{x^4(t)}$$', 'Interpreter', 'Latex');
leg.FontSize=15; leg.Orientation='Horizontal';

xlim([0, 5])
grid on; grid minor;
set(gca, 'Fontsize', 16)
title(['Component 2.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')


sp3 = subplot(2,2,3);
plot(0:5,x1(3:4:end-1), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot(0:5,x2(3:4:end-1), 'Color', '#6593A6', 'Linewidth', 2); 
plot(0:5,x3(3:4:end-1), 'Color', '#A9BF04', 'Linewidth', 2); 
plot(0:5,x4(3:4:end-1), 'Color', '#D95F69', 'Linewidth', 2); 

xlabel('$\textbf{Time}$', 'Interpreter', 'Latex', 'Fontsize', 20);
ylabel(['$$\mathbf{x^i(t)}$$'], 'Interpreter', 'Latex', 'Fontsize', 60);

leg = legend('$$\mathbf{x^1(t)}$$','$$\mathbf{x^2(t)}$$','$$\mathbf{x^3(t)}$$','$$\mathbf{x^4(t)}$$', 'Interpreter', 'Latex');
leg.FontSize=15; leg.Orientation='Horizontal';

xlim([0, 5])
grid on; grid minor;
set(gca, 'Fontsize', 16)
title(['Component 3.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')


sp4 = subplot(2,2,4);
plot(0:5,x1(4:4:end), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot(0:5,x2(4:4:end), 'Color', '#6593A6', 'Linewidth', 2); 
plot(0:5,x3(4:4:end), 'Color', '#A9BF04', 'Linewidth', 2); 
plot(0:5,x4(4:4:end), 'Color', '#D95F69', 'Linewidth', 2); 

xlabel('$\textbf{Time}$', 'Interpreter', 'Latex', 'Fontsize', 20);

leg = legend('$$\mathbf{x^1(t)}$$','$$\mathbf{x^2(t)}$$','$$\mathbf{x^3(t)}$$','$$\mathbf{x^4(t)}$$', 'Interpreter', 'Latex');
leg.FontSize=15; leg.Orientation='Horizontal';

xlim([0, 5])
grid on; grid minor;
set(gca, 'Fontsize', 16)
title(['Component 4.'],...
    'Interpreter', 'Latex', 'Fontsize', 16, 'Fontweight', 'normal')


sgtitle(['\textbf{States Evolution for Projected Subgradient Method.}'],...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold')

sp1.Position = [0.0814 0.5728 0.4239 0.2996];
sp2.Position = [0.5548 0.5728 0.4239 0.2996];
sp3.Position = [0.0814 0.1159 0.4239 0.2996];
sp4.Position = [0.5548 0.1159 0.4239 0.2996];

%% Single Plot
close all
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])

plot(rvec,ivec, 'Color', '#F29544', 'Linewidth', 2); hold on;

xl = xlabel('$$\textbf{Penalty Parameter}$$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{Iterations}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

% leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');
leg = legend('Iterations', '$$\rho = 4$$', '$$\rho = 7$$', '$$\rho = 12$$','Interpreter', 'Latex');

leg.FontSize=30;

tit = title('\textbf{Number of Iterations needed for $$\Delta x_f \leq \varepsilon$$. Varying $$\rho$$.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(rvec)])
grid on; grid minor
set(gca, 'Fontsize', 16)
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.