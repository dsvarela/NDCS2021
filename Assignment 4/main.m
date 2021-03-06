%%
clearvars
close all
clc

%% General Formulation

% From aircraft.m File
Tfinal=5;
umax=100;

% Definition of system dimension
dim.nx=4;     %state dimension
dim.nu=2;     %input dimension
dim.N=Tfinal; %horizon

%Definition of quadratic cost function (same for all aircraft)
weight.Q=eye(dim.nx);   %weight on output
weight.R=eye(dim.nu);   %weight on input

% Anonymous function used to normalize vectors. Will come in useful at some
% point.
nv = @(v) v/norm(v);
%% Problem Formulation for Aircraft 1

% From aircraft.m File
LTI.A=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
LTI.B=[2 0;0 2;3 0;0 3];
LTI.x0=[-10;10;-1;1];

% Generation of prediction model and respective cost
% x = Tx0 + Su
[T1,S1]=predmodgen(LTI,dim);
% J = x'x + u'u = (Tx0 + Su)'(Tx0 + Su) + u'u
% J = (2S'Tx0)u + u'(S'S + I)u = u'Hu + hu
[H1,h1]=costgen(T1,S1,dim,LTI);

% Structure: Aeq1*u + beq1*x0 == xf1
beq1 = T1(dim.nx*(dim.N-1)+1:end,:)*LTI.x0;
Aeq1 = S1(dim.nx*(dim.N-1)+1:end,:);

%% Problem Formulation for Aircraft 2

% From aircraft.m File
LTI.A=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
LTI.B=[3 0; 0 3; 7 0; 0 7];
LTI.x0=[10;10;1;1];

% Generation of prediction model and respective cost
% x = Tx0 + Su
[T2,S2]=predmodgen(LTI,dim);
% J = x'x + u'u = (Tx0 + Su)'(Tx0 + Su) + u'u
% J = (2S'Tx0)u + u'(S'S + I)u = u'Hu + hu
[H2,h2]=costgen(T2,S2,dim,LTI);

% Structure: Aeq2*u + beq2*x0 == xf2
beq2 = T2(dim.nx*(dim.N-1)+1:end,:)*LTI.x0;
Aeq2 = S2(dim.nx*(dim.N-1)+1:end,:);

%% Problem Formulation for Aircraft 3
% From aircraft.m File
LTI.A=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
LTI.B=[1 0; 0 1; 1.1 0; 0 1.1];
LTI.x0=[10;-10;1;-1];

% Generation of prediction model and respective cost
% x = Tx0 + Su
[T3,S3]=predmodgen(LTI,dim);
% J = x'x + u'u = (Tx0 + Su)'(Tx0 + Su) + u'u
% J = (2S'Tx0)u + u'(S'S + I)u = u'Hu + hu
[H3,h3]=costgen(T3,S3,dim,LTI);

% Structure: Aeq3*u + beq3*x0 == xf3
beq3 = T3(dim.nx*(dim.N-1)+1:end,:)*LTI.x0;
Aeq3 = S3(dim.nx*(dim.N-1)+1:end,:);

%% Problem Formulation for Aircraft 4

% From aircraft.m File
LTI.A=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
LTI.B=[6 0;0 6;20 0; 0 20];
LTI.x0=[-10;-10;-1;-1];

% Generation of prediction model and respective cost
% x = Tx0 + Su
[T4,S4]=predmodgen(LTI,dim);
% J = x'x + u'u = (Tx0 + Su)'(Tx0 + Su) + u'u
% J = (2S'Tx0)u + u'(S'S + I)u = u'Hu + hu
[H4,h4]=costgen(T4,S4,dim,LTI);

% Structure: Aeq4*u + beq4*x0 == xf4
beq4 = T4(dim.nx*(dim.N-1)+1:end,:)*LTI.x0;
Aeq4 = S4(dim.nx*(dim.N-1)+1:end,:);

%% General Inequality Constraints

% |ui| <= umax/Tfinal => ui <= umax/Tfinal & -ui <= umax/Tfinal
Ast = [-1;1]; Ain = kron(eye(dim.N*dim.nu), Ast);
bst = [umax/Tfinal; umax/Tfinal]; bin = repmat(bst,dim.N*dim.nu,1);


%% Centralized solution

Hc = blkdiag(H1,H2,H3,H4);
hc = [h1, h2, h3, h4];
Aeqc = [Aeq1, -Aeq2, zeros(dim.nx,dim.nu*dim.N), zeros(dim.nx,dim.nu*dim.N);
    zeros(dim.nx,dim.nu*dim.N), Aeq2, -Aeq3, zeros(dim.nx,dim.nu*dim.N);
    zeros(dim.nx,dim.nu*dim.N), zeros(dim.nx,dim.nu*dim.N), Aeq3, -Aeq4];
beqc = [beq1-beq2;beq2-beq3;beq3-beq4];
Ainc = [blkdiag(Ain,Ain,Ain,Ain)];
binc = ones(2*4*dim.nu*dim.N,1)*umax/Tfinal;

[uc,fval,eflag,~,lambdac] = quadprog(2*Hc,hc',Ainc,binc,Aeqc,-beqc);
x1 = T1*[-10;10;-1;1] + S1*uc(1:10); x1 = [[-10;10;-1;1]; x1];
x2 = T2*[10;10;1;1] + S2*uc(11:20); x2 = [[10;10;1;1]; x2];
x3 = T3*[10;-10;1;-1] + S3*uc(21:30); x3 = [[10;-10;1;-1]; x3];
x4 = T4*[-10;-10;-1;-1] + S4*uc(31:40); x4 = [[-10;-10;-1;-1]; x4];
xfc = x1(end-3:end);

%% Individual Optimization Problems
%initialize
% Constraints on xfi = xfj
l12 = zeros(4,1); l23 = zeros(4,1); l34 = zeros(4,1); l41 = zeros(4,1);

% Constraints on ui <= umax/Tfinal and -ui <= umax/Tfinal
mu1 = zeros(20,1); mu2 = zeros(20,1); mu3 = zeros(20,1); mu4 = zeros(20,1);

xf1v = zeros(4,1); xf2v = zeros(4,1); xf3v = zeros(4,1); xf4v = zeros(4,1);

xavg = ([-10;10;-1;1] + [10;10;1;1] + [10;-10;1;-1] + [-10;-10;-1;-1])/3;
xf1 = xavg; xf2 = xavg; xf3 = xavg; xf4 = xavg;

a = 1;
res = 1e-3;
Deltaxf = Inf;
i = 1;
while Deltaxf > res
    opts = optimoptions('quadprog', 'Display', 'off');
    %Node 1
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    h1d = h1 + (l12'-l41')*Aeq1 + mu1'*Ain;
    % u1 = quadprog(2*H1,h1d,[],[],[],[],[],[],[],opts);
    u1 = -H1\h1d'/2;
    % Store variables for plotting and stopping purposes
    xf1 = Aeq1*u1+beq1; xf1 = xf1(:); xf1v(:,i) = xf1;
    u1v(:,i) = u1;
    % Update l41 with info from node 4
    l41 = l41 + a*(xf4-xf1);
    % Update mu1 with own info
    mu1 = max(mu1 + a*(Ain*u1-bin),0);
    %Node 2
    % Solve Subproblem with info from neighbouring nodes, 1 and 3
    h2d = h2 + (l23'-l12')*Aeq2 + mu2'*Ain;
    % u2 = quadprog(2*H2,h2d,[],[],[],[],[],[],[],opts);
    u2 = -H2\h2d'/2;
    % Store variables for plotting and stopping purposes
    xf2 = Aeq2*u2+beq2; xf2 = xf2(:); xf2v(:,i) = xf2;
    u2v(:,i) = u2;
    % Update l12 with info from node 1
    l12 = l12 + a*(xf1-xf2); %
    % Update mu2 with own info
    mu2 = max(mu2 + a*(Ain*u2-bin),0);
    %Node 3
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    h3d = h3 + (l34'-l23')*Aeq3 + mu3'*Ain;
    % u3 = quadprog(2*H3,h3d,[],[],[],[],[],[],[],opts);
    u3 = -H3\h3d'/2;
    % Store variables for plotting and stopping purposes
    xf3 = Aeq3*u3+beq3; xf3 = xf3(:); xf3v(:,i) = xf3;
    u3v(:,i) = u3;
    % Update l23 with info from node 2
    l23 = l23 + a*(xf2-xf3);
    % Update mu3 with own info
    mu3 = max(mu3 + a*(Ain*u3-bin),0);
    %Node 4
    % Solve Subproblem with info from neighbouring nodes, 3 and 1
    h4d = h4 + (l41'-l34')*Aeq4 + mu4'*Ain;
    % u4 = quadprog(2*H4,h4d,[],[],[],[],[],[],[],opts);
    u4 = -H4\h4d'/2;
    % Store variables for plotting and stopping purposes
    xf4 = Aeq4*u4+beq4; xf4 = xf4(:); xf4v(:,i) = xf4;
    u4v(:,i) = u4;
    % Update l34 with info from node 3
    l34 = l34 + a*(xf3-xf4);
    % Update mu4 with own info
    mu4 = max(mu4 + a*(Ain*u4-bin),0);
    
    % Store dual variables for plotting purposes
    ls(:,i) = [l12;l23;l34;l41];
    mus(:,i) = [mu1;mu2;mu3;mu4];
    
    % Check for convergence.
    Deltaxf = max([norm(xf1-xfc), norm(xf2-xfc), norm(xf3-xfc), norm(xf4-xfc)])
    
    i = i+1
end

%% Plotting Results
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])

plot((vecnorm(xf1v-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2);
plot((vecnorm(xf3v-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2);
plot((vecnorm(xf4v-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2);

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||x_f^i - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');

leg.FontSize=30;

tit = title('\textbf{Projected Subgradient Method - Error Evolution Normalized Against Centralized Solution.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xf1v)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.

%% Individual Optimization Problems (Nesterov's Method)

Aeq_ = [Aeq1 -Aeq2 zeros(dim.nx,dim.nu*dim.N) zeros(dim.nx,dim.nu*dim.N);
    zeros(dim.nx,dim.nu*dim.N) Aeq2 -Aeq3 zeros(dim.nx,dim.nu*dim.N);
    zeros(dim.nx,dim.nu*dim.N) zeros(dim.nx,dim.nu*dim.N) Aeq3 -Aeq4;
    -Aeq1 zeros(dim.nx,dim.nu*dim.N) zeros(dim.nx,dim.nu*dim.N) Aeq4];
beq_ = [beq2 - beq1; beq3 - beq2; beq4 - beq3; beq1 - beq4];

Ain_ = blkdiag(Ain, Ain, Ain, Ain);
bin_ = repmat(bin, 4,1);

Anest = [Aeq_; Ain_]; bnest = [beq_; bin_];
Hnest = blkdiag(H1,H2,H3,H4); hnest = [h1, h2, h3, h4]';

L = norm(Anest/Hnest*Anest'/2,2);
%initialize
% Constraints on xfi = xfj
l12 = zeros(4,1); l23 = zeros(4,1); l34 = zeros(4,1); l41 = zeros(4,1);
ls = [[l12;l23;l34;l41] , [l12;l23;l34;l41]];
% Constraints on ui <= umax/Tfinal and -ui <= umax/Tfinal
mu1 = zeros(20,1); mu2 = zeros(20,1); mu3 = zeros(20,1); mu4 = zeros(20,1);
mus = [[mu1;mu2;mu3;mu4] , [mu1;mu2;mu3;mu4]];


xf1v = zeros(4,1); xf2v = zeros(4,1); xf3v = zeros(4,1); xf4v = zeros(4,1);

xavg = ([-10;10;-1;1] + [10;10;1;1] + [10;-10;1;-1] + [-10;-10;-1;-1])/3
xf1 = xavg; xf2 = xavg; xf3 = xavg; xf4 = xavg;

Deltaxf = inf;
i= 1;
while Deltaxf > res
    opts = optimoptions('quadprog', 'Display', 'off');
    %Problem 1
    h1d = h1 + (l12'-l41')*Aeq1 + mu1'*Ain;
    u1 = quadprog(2*H1,h1d,[],[],[],[],[],[],[],opts);
    xf1 = Aeq1*u1+beq1; xf1 = xf1(:); xf1v(:,i) = xf1;
    %Problem 2
    h2d = h2 + (l23'-l12')*Aeq2 + mu2'*Ain;
    u2 = quadprog(2*H2,h2d,[],[],[],[],[],[],[],opts);
    xf2 = Aeq2*u2+beq2; xf2 = xf2(:); xf2v(:,i) = xf2;
    %Problem 3
    h3d = h3 + (l34'-l23')*Aeq3 + mu3'*Ain;
    u3 = quadprog(2*H3,h3d,[],[],[],[],[],[],[],opts);
    xf3 = Aeq3*u3+beq3; xf3 = xf3(:); xf3v(:,i) = xf3;
    %Problem 4
    h4d = h4 + (l41'-l34')*Aeq4 + mu4'*Ain;
    u4 = quadprog(2*H4,h4d,[],[],[],[],[],[],[],opts);
    xf4 = Aeq4*u4+beq4; xf4 = xf4(:); xf4v(:,i) = xf4;
    %Dual variables update
    zm1 = [ls(:,i); mus(:,i)];
    z = [ls(:,i+1); mus(:,i+1)];
    
    
    % Tuning parameter
    beta = 0.75;
    
    % 	 alph = (1-sqrt(i))/(1+sqrt(i));
    alph = beta*(i-1)/(i+2);
    
    vk = z + alph*(z-zm1);
    
    gradf = 1/2 * Anest/Hnest*(Anest'*vk + hnest) + bnest;
    
    zp1 = vk - 1/L*gradf;
    
    ls(:,i+2) = zp1(1:16); % ls = [l-1 , l0, l1, l2, ...]
    mus(:,i+2) = max(zp1(17:end), 0); % mus = [mu-1, mu0, mu1, mu2, ...]
    
    l12 = ls(1:4,i+2); l23 = ls(5:8,i+2); l34 = ls(9:12,i+2); l41 = ls(13:16,i+2);
    mu1 = mus(1:20,i+2); mu2 = mus(21:40,i+2); mu3 = mus(41:60,i+2); mu4 = mus(61:80,i+2);
    
    % Check for convergence.
    Deltaxf = max([norm(xf1-xf2), norm(xf2-xf3), norm(xf3-xf4), norm(xf4-xf1)])
    i = i+1
    
    % a = a/i;
    % If it's going to blow up
    if Deltaxf > 10
        break
    end
    
end

%% Plotting Results
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])

plot((vecnorm(xf1v-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2);
plot((vecnorm(xf3v-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2);
plot((vecnorm(xf4v-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2);

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||x_f^i - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');

leg.FontSize=30;

tit = title('\textbf{Nesterov Method - Error Evolution Normalized Against Centralized Solution.}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xf1v)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.

%% Consensus Iterations
xf1 = zeros(4,1); xf2 = zeros(4,1); xf3 = zeros(4,1); xf4 = zeros(4,1);
xf1v = zeros(4,1); xf2v = zeros(4,1); xf3v = zeros(4,1); xf4v = zeros(4,1);
% xf1 = xfc; xf2 = xfc; xf3 = xfc; xf4 = xfc;

W = [0.75 0.25 0.00 0.00;
    0.25 0.50 0.25 0.00;
    0.00 0.25 0.50 0.25;
    0.00 0.00 0.25 0.75;];
Wavg = 1/4*ones(4,4);
Deltaxf = Inf;

a0 = 0.05;
a = a0;
i = 2;
res = 1e-4;
while Deltaxf > res
    opts = optimoptions('quadprog', 'Display', 'off');
    %Node 1
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    [u1,~,~,~,l1] = quadprog(2*H1,h1,Ain,bin,Aeq1,-beq1+xf1,[],[],[],opts);
    % Store variables for plotting and stopping purposes
    u1v(:,i) = u1; xf1v(:,i) = xf1;
    xf1loc = xf1 + a*l1.eqlin;
    %Node 2
    % Solve Subproblem with info from neighbouring nodes, 1 and 3
    [u2,~,~,~,l2] = quadprog(2*H2,h2,Ain,bin,Aeq2,-beq2+xf2,[],[],[],opts);
    % Store variables for plotting and stopping purposes
    u2v(:,i) = u2; xf2v(:,i) = xf2;
    xf2loc = xf2 + a*l2.eqlin;
    %Node 3
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    [u3,~,~,~,l3] = quadprog(2*H3,h3,Ain,bin,Aeq3,-beq3+xf3,[],[],[],opts);
    % Store variables for plotting and stopping purposes
    u3v(:,i) = u3; xf3v(:,i) = xf3;
    xf3loc = xf3 + a*l3.eqlin;
    %Node 4
    % Solve Subproblem with info from neighbouring nodes, 3 and 1
    [u4,~,~,~,l4] = quadprog(2*H4,h4,Ain,bin,Aeq4,-beq4+xf4,[],[],[],opts);
    % Store variables for plotting and stopping purposes
    u4v(:,i) = u4; xf4v(:,i) = xf4;
    xf4loc = xf4 + a*l4.eqlin;
    
    sub = [xf1loc'; xf2loc'; xf3loc'; xf4loc'];
    
    update = W^20*sub;
    xf1 = update(1,:)';
    xf2 = update(2,:)';
    xf3 = update(3,:)';
    xf4 = update(4,:)';
    
    Deltaxf = max([norm(xf1-xf1v(:,end-1)),...
        norm(xf2-xf2v(:,end-1)), norm(xf3-xf3v(:,end-1)), norm(xf4-xf4v(:,end-1))]);
    
    i = i+1
    % a = a0/(i-1);
    % No convergence if the alpha is varying, so break after 500
    % iterations.
    if i == 500
        break
    end
end
%% Plotting Results
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.725]);
% figure('Position', [10, 50, 1400, 600])

plot((vecnorm(xf1v-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;
plot((vecnorm(xf2v-xfc)/norm(xfc)), 'Color', '#6593A6', 'Linewidth', 2);
plot((vecnorm(xf3v-xfc)/norm(xfc)), 'Color', '#A9BF04', 'Linewidth', 2);
plot((vecnorm(xf4v-xfc)/norm(xfc)), 'Color', '#D95F69', 'Linewidth', 2);

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||x_f^i - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

leg = legend('$$\mathbf{x_f^1}$$','$$\mathbf{x_f^2}$$','$$\mathbf{x_f^3}$$','$$\mathbf{x_f^4}$$', 'Interpreter', 'Latex');

leg.FontSize=30;

tit = title({'\textbf{Combined Consensus/Incremental Subgradient}','\textbf{Error Evolution Normalized Against Centralized Solution.}','$$\alpha = 0.05,\quad \varphi = 20$$'},...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xf1v)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.
%% Dual Ascent
% xf1 = zeros(4,1); xf2 = zeros(4,1); xf3 = zeros(4,1); xf4 = zeros(4,1);
% xf1v = zeros(4,1); xf2v = zeros(4,1); xf3v = zeros(4,1); xf4v = zeros(4,1);
% % xf1 = xfc; xf2 = xfc; xf3 = xfc; xf4 = xfc;
%
% l1 = zeros(4,1); l2 = zeros(4,1); l3 = zeros(4,1); l4 = zeros(4,1);
% ls = zeros(16,1);
%
% Deltaxf = Inf; Deltal = Inf;
%
% a = 0.1; al= 1;
% i = 2;
% while Deltaxf > -1e-4
%     opts = optimoptions('quadprog', 'Display', 'off');
%     j=1;
% while Deltal > -1
%     %Node 1
%     h1d = h1 + l1'*Aeq1;
%     [u1,~,~,~,l1f] = quadprog(2*H1,h1d,Ain,bin,[],[],[],[],[],opts);
%     %Node 2
%     h2d = h2 + l2'*Aeq2;
%     [u2,~,~,~,l2f] = quadprog(2*H2,h2d,Ain,bin,[],[],[],[],[],opts);
%     %Node 3
%     % Solve Subproblem with info from neighbouring nodes, 2 and 4
%     h3d = h3 + l3'*Aeq3;
%     [u3,~,~,~,l3f] = quadprog(2*H3,h3d,Ain,bin,[],[],[],[],[],opts);
%     %Node 4
%     h4d = h4 + l4'*Aeq4;
%     [u4,~,~,~,l4f] = quadprog(2*H4,h4d,Ain,bin,[],[],[],[],[],opts);
%
%     l1 = l1+al*(Aeq1*u1 + beq1 - xf1);
%     l2 = l2+al*(Aeq2*u2 + beq2 - xf2);
%     l3 = l3+al*(Aeq3*u3 + beq3 - xf3);
%     l4 = l4+al*(Aeq4*u4 + beq4 - xf4)
%     ls(:,j) = [l1;l2;l3;l4];
%     j = j+1;
%     if j == 1000
%         j;
%     end
% end
%     % Store variables for plotting and stopping purposes
%     u1v(:,i) = u1; xf1v(:,i) = xf1;
%     u2v(:,i) = u2; xf2v(:,i) = xf2;
%     u3v(:,i) = u3; xf3v(:,i) = xf3;
%     u4v(:,i) = u4; xf4v(:,i) = xf4;
%     xf1loc = xf1 + a*l1;
%     xf2loc = xf2 + a*l2;
%     xf3loc = xf3 + a*l3;
%     xf4loc = xf4 + a*l4;
%
%     sub = [xf1loc'; xf2loc'; xf3loc'; xf4loc'];
%     update = W^100*sub;
%     xf1 = update(1,:)';
%     xf2 = update(2,:)';
%     xf3 = update(3,:)';
%     xf4 = update(4,:)';
%
%     Deltaxf = max([norm(xf1-xf1v(:,end-1)),...
%         norm(xf2-xf2v(:,end-1)), norm(xf3-xf3v(:,end-1)), norm(xf4-xf4v(:,end-1))]);
%     i = i+1;
%     a = 0.5/i;
%     % If it's going to blow up
%     if i == 1000
%         break
%     end
%
% end

%% ADMM + Consensus Iterations


xf = zeros(4,1); xfv = zeros(4,1);

l1 = zeros(4,1); l2 = zeros(4,1); l3 = zeros(4,1); l4 = zeros(4,1);
ls = zeros(16,1);


rho = 6.8;
i = 1;
Deltaxf = Inf;
while Deltaxf > 1e-4
    opts = optimoptions('quadprog', 'Display', 'off');
    %Node 1
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    H1d = (H1 + rho/2*Aeq1'*Aeq1) + (H1 + rho/2*Aeq1'*Aeq1)';
    u1 = quadprog(H1d,(h1 + l1'*Aeq1 + rho*(beq1 - xf)'*Aeq1),Ain,bin,[],[],[],[],[],opts);
    % Store variables for plotting and stopping purposes
    
    %Node 2
    % Solve Subproblem with info from neighbouring nodes, 1 and 3
    H2d = (H2 + rho/2*Aeq2'*Aeq2)+  (H2 + rho/2*Aeq2'*Aeq2)';
    u2 = quadprog(H2d,(h2 + l2'*Aeq2 + rho*(beq2 - xf)'*Aeq2),Ain,bin,[],[],[],[],[],opts);
    % Store variables for plotting and stopping purposes
    
    %Node 3
    % Solve Subproblem with info from neighbouring nodes, 2 and 4
    H3d = (H3 + rho/2*Aeq3'*Aeq3) + (H3 + rho/2*Aeq3'*Aeq3)';
    u3 = quadprog(H3d,(h3 + l3'*Aeq3 + rho*(beq3 - xf)'*Aeq3),Ain,bin,[],[],[],[],[],opts);
    % Store variables for plotting and stopping purposes
    %Node 4
    % Solve Subproblem with info from neighbouring nodes, 3 and 1
    H4d = (H4 + rho/2*Aeq4'*Aeq4) + (H4 + rho/2*Aeq4'*Aeq4)';
    u4 = quadprog(H4d,(h4 + l4'*Aeq4 + rho*(beq4 - xf)'*Aeq4),Ain,bin,[],[],[],[],[],opts);
    % Store variables for plotting and stopping purposes
    
    xf = ((l1/rho + Aeq1*u1 + beq1) + (l2/rho + Aeq2*u2 + beq2)+ ...
        (l3/rho + Aeq3*u3 + beq3) + (l4/rho + Aeq4*u4 + beq4))/4;
    
    l1 = l1+rho*(Aeq1*u1 + beq1 - xf);
    l2 = l2+rho*(Aeq2*u2 + beq2 - xf);
    l3 = l3+rho*(Aeq3*u3 + beq3 - xf);
    l4 = l4+rho*(Aeq4*u4 + beq4 - xf);
    ls(:,i) = [l1;l2;l3;l4];
    
    xfv(:,i+1) = xf;
    Deltaxf = norm(xf-xfv(:,end-1));
    
    i = i+1
    
    if i == 1000
        fprintf('Here.')
        break
    end
end

%% Plotting Results
figure('Position', [10, 50, 1400, 600],'DefaultAxesPosition', [0.085, 0.125, 0.875, 0.775]);
% figure('Position', [10, 50, 1400, 600])

plot((vecnorm(xfv-xfc)/norm(xfc)), 'Color', '#F29544', 'Linewidth', 2); hold on;

xl = xlabel('$\textbf{Iterations}$', 'Interpreter', 'Latex', 'Fontsize', 20);
xl.FontSize=50;

yl = ylabel('$$\mathbf{\frac{||x_f - x_f^*||_2}{||x_f^*||}}$$', 'Interpreter', 'Latex', 'Fontsize', 50);

leg = legend('Consensus Variable $$\mathbf{x_f}$$', 'Interpreter', 'Latex');
leg.FontSize=30;

tit = title('\textbf{ADMM - Error Evolution Normalized Against Centralized Solution. $$\rho = 6.8$$}',...
    'Interpreter', 'Latex', 'Fontsize', 20, 'Fontweight', 'bold');
tit.FontSize=20;

xlim([0,length(xfv)])
grid on; grid minor
set(gca, 'Fontsize', 16, 'Yscale', 'log')
ax = gca;
ax.GridAlpha = 0.25;  % Make grid lines less transparent.
