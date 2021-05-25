clearvars; clc; close all;
%% Student Number
SN = [5 3 2 1 2 6 3];
a = SN(1); b = SN(3); c = SN(end);
syms h s t hl
%% Matrices
A = [a b+0.5; 0 -c];
B = [0; 1];

%% Find Static Controller
P = [-2,-3];
K = place(A,B,P); % Find controller with the Place Command
AK = A-B*K; % Close the loop

%% Find matrices
Fx = expm(A*h); % Fx = e^(Ah)
G = int(expm(A*s)*[0;1],s,0,h); % G = int_{0}^{h} e^(As) ds

FxK = Fx-G*K; % Closed loop

Fxf = matlabFunction(Fx);
Gf = matlabFunction(G);

Fcl = expm(A*hl) - expm(A*(hl-h))*int(expm(A*s)*B*K,s,0,h);
Fclf = matlabFunction(Fcl);

%% Question 3.1 - Max consecutive dropouts with to-zero mechanism
% In effect, this will be an evergrowing switched system, with dynamics:
% (F(h)-G(h)*K)*(F(h)^del)
% with del the number of consecutive dropouts, del = 0,1,...,N

hvec = 0:0.005:0.35;
del_ = zeros(1,length(hvec));

for j = 1:length(hvec)
    h = hvec(j);
    Fh = Fxf(h);
    Gh = Gf(h); % Store correponding matrices
    CL = Fh - Gh*K; % close the loop
    
    % Initialize Things
    del = 0;
    % hl = [h, 2*h, 3*h, 4*h, 5*h, 6*h];
    n = 2; % Matrix size
    while 1
        cvx_begin sdp
        variable M(n,n) symmetric
        for i = 0:del
            % [-M, M*Fclf(h,hl(i+1))'; Fclf(h,hl(i+1))*M, -M] <= -eye(2*n);
            [-M, M*(CL*(Fh^i))'; (CL*(Fh^i))*M, -M] <= -eye(2*n);
        end
        cvx_end
        % P = inv(M);
        if strcmp(cvx_status,'Solved')
            del = del+1;
        else
            del_(j) = max(del-1,0);
            break
        end
    end
end
% Plot Things
figure(1);
set(gcf,'Position',[10 50 1210 450]);
scatter(hvec, del_ , [],  [238, 0, 0]/256, 'o', 'filled'); hold on;
plot(hvec, del_,'k--', 'Linewidth',0.1); hold off;
legend('max \delta', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Maximum Number of Consecutive Packet Dropouts \delta with Sampling Time {\it\fontname{Cambria Math}h}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\delta', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;

%% 3.1 Simulation Test Bench
% Brute force a couple hundred simulations, see what comes out.

h = 0.1;
Fh = Fxf(h);
Gh = Gf(h); % Store correponding matrices
clear CL
simlen = 100
x1 = zeros(250,simlen+1);
x2 = zeros(250,simlen+1);
for j = 1:250
    x0 = rand(2,1)*20-10;
    x0(1) = abs(x0(1))+10; x0(2) = -abs(x0(2))-10;
    x = zeros(2,simlen+1); x(:,1) = x0;
    x1(j,1) = x0(1);
    x2(j,1) = x0(2);
    
    for i = 1:simlen
        if rem(i,9) == 0 || rem(i,9) == 1 || rem(i,9) == 2
            A_ = Fh;
        else
            A_ = Fh - Gh*K;
        end
        x(:,i+1) = A_*x(:,i);
        x1(j,i+1) = x(1,i+1);
        x2(j,i+1) = x(2,i+1);
    end
end

x1avg = mean(x1);
x2avg = mean(x2);

figure(2);
set(gcf,'Position',[10 50 1210 450]);
plot((0:simlen)*h, x1avg,'Color', '#EE0000', 'Linewidth',1.5);  hold on;
plot((0:simlen)*h, x2avg,'Color', '#EE0000', 'Linewidth',1.5);

plot((0:simlen)*h, x1,'Color', '#BBBBBB', 'Linewidth',0.1);
plot((0:simlen)*h, x2,'Color', '#BBBBBB', 'Linewidth',0.1);
plot((0:simlen)*h, x1avg,'Color', '#EE0000', 'Linewidth',1.5);
plot((0:simlen)*h, x2avg,'Color', '#EE0000', 'Linewidth',1.5);  hold off;

legend('x̅_1', 'x̅_2', 'x_1', 'x_2', 'Fontweight', 'bold', 'Orientation', 'Horizontal',...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Evolution of states with 1 out of every 3 of packets lost', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('Amplitude', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('Time (s)', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;



%% Question 3.2 - upper bound p* on the probability of a Bernoulli process
% From the lectures, we know MSS is guaranteed if:
% P - (1-p)*A0'*P*A0 - p*A1'*P*A1 > 0
% A0 is the dynamics with no packet loss: A0 = (F(h) - G(h)*K)
% A1 is the dynamics with packet losses:  A1 = F(h)
n = 2;
hvec = 0.01:0.025:0.35;
p_ = 0:0.01:0.3;
for j = 1:length(hvec)
    h = hvec(j)
    Fh = Fxf(h);
    Gh = Gf(h);
    CL = Fh - Gh*K; % close the loop
    for i = 1:length(p_)
        p = p_(i)
        cvx_begin sdp quiet
        variable P(n,n) symmetric
        P - (1-p)*CL'*P*CL - p*Fh'*P*Fh >= eye(n);
        P >= eye(n);
        cvx_end
        if strcmp(cvx_status,'Solved') & p < 1
            pmat(i,j) = 1; % Flag for solved.
        end
    end
end
[rows,cols] = find(pmat == 1);
%  Plot Things
c_ = 0;
j =1;
for i = length(cols):-1:1
    if cols(i) ~= c_
        c_ = cols(i);
        mcols(j) = cols(i);
        mrows(j) = rows(i);
        j = j+1;
    end
end
%%
figure(2);
set(gcf,'Position',[10 50 1210 450]);
scatter(cols*0.025, (rows-1)*0.01, [], [150 150 150]/256, 'x','LineWidth',2); hold on;
scatter(mcols*0.025, (mrows-1)*0.01, [], [0 0 0], 'rx','LineWidth',2)
legend('{\it\fontname{Cambria Math}p}', '{\it\fontname{Cambria Math}p*}', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Evolution of upper bound on {\it\fontname{Cambria Math}p} with Sampling Time {\it\fontname{Cambria Math}h}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('{\it\fontname{Cambria Math}p}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;
hold off
%% Question 3.3
% Brute force a couple hundred simulations, see what comes out.

h = 0.1; pstar = .22; % h = 0.1, and corresponding p*
Fh = Fxf(h);
Gh = Gf(h);
clear CL
CL{1} = Fh - Gh*K; % close the loop A0
CL{2} = Fh; % close the loop A1
simlen = 200;
attempts = 500;
x1 = zeros(attempts,simlen+1);
x2 = zeros(attempts,simlen+1);
xn = zeros(attempts,simlen+1);
for j = 1:attempts
    dyn = randi(2); % Initial dynamics picked at random
    A_ = CL{dyn};
    x0 = rand(2,1)*20-10;
    x = zeros(2,simlen+1); x(:,1) = x0;
    x1(j,1) = x0(1);
    x2(j,1) = x0(2);
    xn(j,1) = norm(x0)^2;
    for i = 1:simlen
        if rand < pstar
            A_ = CL{2};
        else
            A_ = CL{1};
        end
        x(:,i+1) = A_*x(:,i);
        x1(j,i+1) = x(1,i+1);
        x2(j,i+1) = x(2,i+1);
        xn(j,i+1) = norm([x(1,i+1);x(2,i+1)])^2;
    end
end

x1avg = mean(x1);
x2avg = mean(x2);
xnavg = mean(xn);

figure(2);
set(gcf,'Position',[10 50 1210 450]); 
hold on;
plot((0:simlen)*h, xnavg,'Color', '#000000', 'Linewidth',1.5);
plot((0:simlen)*h, x1avg,'Color', '#EE0000', 'Linewidth',1.5);  
plot((0:simlen)*h, x2avg,'Color', '#EE0000', 'Linewidth',1.5);

plot((0:simlen)*h, x1,'Color', '#BBBBBB', 'Linewidth',0.1);
plot((0:simlen)*h, x2,'Color', '#BBBBBB', 'Linewidth',0.1);
plot((0:simlen)*h, x1avg,'Color', '#EE0000', 'Linewidth',1.5);
plot((0:simlen)*h, x2avg,'Color', '#EE0000', 'Linewidth',1.5);  
plot((0:simlen)*h, xnavg,'Color', '#000000', 'Linewidth',1.5);

hold off;
legend('E(||x||^2)', 'Fontweight', 'bold', 'Orientation', 'vertical',...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southeast');

sgtitle('Evolution of states with {\it\fontname{Cambria Math}h}=0.1s for {\it\fontname{Cambria Math}p*} = 0.22', ...
    'Fontsize', 20, 'Fontweight', 'bold')
set(gca, 'Fontsize', 15)
ylabel('E[||x||^2]', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('Time (s)', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;

%% Question 3.4
hvec = 0.1;
p = 0.35;
Fh = Fxf(h);
Gh = Gf(h); % Store correponding matrices
clear CL
CL{1} = Fh - Gh*K; % close the loop A0
CL{2} = Fh; % close the loop A1

n = 2;
lambda = 0.5:0.01:0.6; lambda(1) = []; lambda(end) = [];
for i = 1:length(lambda)
    L = (1/(lambda(i))^(1-p))^(1/p) - eps;
    cvx_begin sdp
    variable P(n,n) symmetric
    % 0.5529
    CL{1}'*P*CL{1} - lambda(i)*P <= zeros(n);
    CL{2}'*P*CL{2} - L*P <= zeros(n);
    P >= eye(n);
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        Lb= L;
        lambdab= lambda(i);
        Pb = P;
        
        break
    end
end

%% Question 3.5 - combinations (p,q) on the probability of a Bernoulli process
% From the lectures, we know MSS is guaranteed if:
% P0 - p*A0'*P0*A0 - (1-p)*A1'*P1*A1 > 0
% P1 - q*A1'*P1*A1 - (1-q)*A0'*P0*A0 > 0
% A0 is the dynamics with no packet loss: A0 = (F(h) - G(h)*K)
% A1 is the dynamics with packet losses:  A1 = F(h)


hvec = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]; % Fix a value, cause I'm tired of doing this for multiple ones.
qmax = [0.625, 0.375, 0.225, 0.125, 0.075,  0.025, 0.025];
pmin = [0.575, 0.525, 0.375, 0.000, 0.000,  0.750, 1];
for k = 1:length(hvec)
    h = hvec(k);
    Fh = Fxf(h);
    Gh = Gf(h); % Store correponding matrices
    CL = Fh - Gh*K; % close the loop
    p_ = pmin(k):0.025:1;
    q_ = 0:0.025:qmax(k);
    for i = 1:length(p_)
        for j = 1:length(q_)
            p = p_(i);
            q = q_(j);
            % Initialize Things
            n = 2; % Matrix size
            cvx_begin sdp quiet
            variable P0(n,n) symmetric semidefinite
            variable P1(n,n) symmetric semidefinite
            P0 - p*CL'*P0*CL - (1-p)*Fh'*P1*Fh >= eye(n);
            P1 - q*Fh'*P1*Fh - (1-q)*CL'*P0*CL >= eye(n);
            cvx_end
            if strcmp(cvx_status,'Solved')
                pq(i,j,k) = 1;
            end
        end
    end
end


% Find rows and columns for which the eigenvalues are smaller than 1
figure(3);
set(gcf,'Position',[10 50 1210 600]);

figure(4);
set(gcf,'Position',[10 50 1210 450]);

for k = 1:length(hvec)
    [rows,cols,~] = find(pq(:,:,k)  == 1);
    p_ = pmin(k):0.025:1;
    q_ = 0:0.025:qmax(k);
    % Plot Things (on the same figure as Question 2.1).
    
    % Also adding back the h to t before plotting, to get the real delay
    % tau = t+h \in [h,2*h)
    figure(3);
    scatter3( ones(1,length(cols))*hvec(k), p_(rows),q_(cols), [], [0 0 0], 'rx','LineWidth',2); hold on;
    % plot3(ones(1,length((0:0.1:1)))*hvec(k), (0:0.1:1),1-(0:0.1:1), 'k--')
    
    figure(4)
    idx = find(p_(rows) == 1 - q_(cols));
    p = p_(rows); p = p(idx);
    scatter(ones(1,length(p))*hvec(k), 1-p, [], [0 0 0], 'rx','LineWidth',2); hold on
end
figure(3);
grid on; grid minor;
legend('{\fontname{Cambria Math}Stable Combinations}',...
    'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');

set(gca, 'Fontsize', 15)
title('Stable Combinations of {\it\fontname{Cambria Math}p_{00}} and {\it\fontname{Cambria Math}p_{11}} in the Gilbert-Elliot Process.', ...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
zlabel('{\it\fontname{Cambria Math}p_{11}}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('{\it\fontname{Cambria Math}p_{00}}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')

figure(4);
grid on; grid minor;
legend('{\fontname{Cambria Math}Stable Combinations}',...
    'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');

set(gca, 'Fontsize', 15)
title('Stable Combinations of where Gilbert-Elliot and Bernoulli are Equivalent.', ...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('{\it\fontname{Cambria Math}p}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
%% Question 3.6

hvec = 0.025:0.025:0.35; % Fix a value, cause I'm tired of doing this for multiple ones.
q_ = 0:0.01:0.3;
pq = zeros(numel(q_),numel(hvec));
for k = 1:length(hvec)
    h = hvec(k);
    Fh = Fxf(h);
    Gh = Gf(h); % Store correponding matrices
    CL = Fh - Gh*K; % close the loop
    for i = 1:length(q_)
        q = q_(i);
        p = 1-q_(i);
        % Initialize Things
        n = 2; % Matrix size
        cvx_begin sdp quiet
        variable P0(n,n) symmetric semidefinite
        variable P1(n,n) symmetric semidefinite
        P0 - p*CL'*P0*CL - (1-p)*Fh'*P1*Fh >= 10^(-3)*eye(n);
        P1 - q*Fh'*P1*Fh - (1-q)*CL'*P0*CL >= 10^(-3)*eye(n);
        P1 >= 10^(-3)*eye(n); P0 >= 10^(-3)*eye(n);
        cvx_end
        if strcmp(cvx_status,'Solved')
            pq(i,k) = 1;
        end
    end
end

[rows,cols] = find(pq == 1);
% Plot Things
c_ = 0;
j =1;
for i = length(cols):-1:1
    if cols(i) ~= c_
        c_ = cols(i);
        mcols(j) = cols(i);
        mrows(j) = rows(i);
        j = j+1;
    end
end

figure(2);
set(gcf,'Position',[10 50 1210 450]);
scatter(cols*0.025, (rows-1)*0.01, [], [150 150 150]/256, 'x','LineWidth',2); hold on;
scatter(mcols*0.025, (mrows-1)*0.01, [], [0 0 0], 'rx','LineWidth',2)
legend('{\it\fontname{Cambria Math}p_{11}}', '{\it\fontname{Cambria Math}p_{11}*}', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Evolution of upper bound on {\it\fontname{Cambria Math}p_{11} = 1-p_{00}} with Sampling Time {\it\fontname{Cambria Math}h}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('{\it\fontname{Cambria Math}p_{11}}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;

%% Question 4.1
syms h s t hl
M0  = int(expm(A*s)*[0;1],s,0,h-t); % M0 = int_{0}^{h-tau} e^(As) ds (G1)
M1  = int(expm(A*s)*[0;1],s,h-t,h); % M1 = int_{h-tau}^{h} e^(As) ds (Fu)
Fe = simplify([[Fx, M1]; zeros(1,3)]);
Ge = simplify([M0; 1]);

alpha1 = exp(-3*(h-t)); alpha1f = matlabFunction(alpha1);
alpha2 = exp(5*(h-t)); alpha2f = matlabFunction(alpha2);

F0 = matlabFunction([exp(5*h), 5/16*(exp(5*h)-exp(-3*h)), 5/48*exp(-3*h)+1/16*exp(5*h);...
    0        , exp(-3*h)                , -1/3*(exp(-3*h));...
    0        , 0                        , 0]);
F1 = [0, 0, -5/48; 0, 0, 1/3; 0, 0, 0];
F2 = [0, 0, -1/16; 0, 0, 0;   0, 0, 0];
% Ftil = simplify(F0 + alpha1*F1 + alpha2*F2);

G0 = [-1/6; 1/3; 1];
G1 = [5/48; -1/3; 0];
G2 = [1/16; 0; 0];
% Gtil = simplify(G0 + alpha1*G1 + alpha2*G2);

%%
hvec= 0.01:0.01:0.10;
mvec = [linspace(0.02,0.105,10), linspace(0.11,0.01,37)];
[rows, cols] = find(goodpts==1);
for j = 1:numel(hvec)
idx = find(rows == j);
colvec = cols(idx);
h = hvec(j);
tauM = (max(colvec)-1)*0.001;
taum = (min(colvec)-1)*0.001;
% taum = 0; tauM = 0.1;
tvec = linspace(taum,tauM,100);

clear a1vec a2vec
for i=1:length(tvec)
    a1vec(i) = alpha1f(h,tvec(i));
    a2vec(i) = alpha2f(h,tvec(i));
end


alphasq= [alpha1f(h,tauM), alpha2f(h,tauM);...
          alpha1f(h,taum), alpha2f(h,tauM);...
          alpha1f(h,taum), alpha2f(h,taum);...
          alpha1f(h,tauM), alpha2f(h,taum)];
m = mvec(j);
f = 10;
alphahex = [alpha1f(h,taum), alpha2f(h,taum);
              (alpha1f(h,taum)*(1-m/f)+alpha1f(h,tauM)*m/f), alpha2f(h,taum);...
              alpha1f(h,tauM), (alpha2f(h,taum)*m/f+alpha2f(h,tauM)*(1-m/f));...
              alpha1f(h,tauM), alpha2f(h,tauM);...
              (alpha1f(h,taum)*m+alpha1f(h,tauM)*(1-m)), alpha2f(h,tauM);...
              alpha1f(h,taum), (alpha2f(h,taum)*(1-m)+alpha2f(h,tauM)*m)];
figure(4);
set(gcf,'Position',[10 50 1210 450]);

plot(a1vec,a2vec,'Color','#EE0000','Linewidth',2); hold on;
sq = polyshape(alphasq(:,1),alphasq(:,2));
sqpl = plot(sq); sqpl.FaceAlpha = 0.1; sqpl.EdgeColor= '#00EE00';
sqpl.FaceColor= '#00EE00'; sqpl.LineWidth = 1.5;
%hex = polyshape(alphahex(:,1),alphahex(:,2));
%hexpl = plot(hex); hexpl.EdgeColor= '#0000EE';
%hexpl.FaceColor= '#0000EE'; hexpl.FaceAlpha = 0.1; hexpl.LineWidth = 1.5;
scatter(alphasq(:,1), alphasq(:,2), [], [0,238,0]/256, 'x', 'LineWidth', 2);
% scatter(alphahex(:,1), alphahex(:,2),[], [0,0,238]/256, 'x', 'LineWidth', 2);
plot(a1vec,a2vec,'Color','#EE0000','Linewidth',2); hold on;
grid on; grid minor;
set(gca, 'Fontsize', 15)
title('Evolution of \alpha_1 and \alpha_2 with the delay \tau for a fixed sampling interval {\it\fontname{Cambria Math}h}=0.1s.', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\alpha_2', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold');
xlabel('\alpha_1', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold');
legend('\alpha_1 vs \alpha_2', 'Initial Polygon', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
xlim([alpha1f(h,taum)-0.005,alpha1f(h,tauM)+0.005])
ylim([alpha2f(h,tauM)-0.05,alpha2f(h,taum)+0.05])
hold off;
end
%% Question 4.2
h = 0.3; tauM = 0.028; taum = 0;
alpha1lim = [alpha1f(h,taum) alpha1f(h,tauM)];
alpha2lim= [alpha2f(h,taum) alpha2f(h,tauM)];


K = [22.4 7 0];
clear CL
Ftil{length(alpha1lim)+length(alpha2lim)} = [];
Gtil{length(alpha1lim)+length(alpha2lim)} = [];
CLtil{length(alpha1lim)+length(alpha2lim)} = [];
k = 1;
for i=1:length(alpha1lim)
    for j = 1:length(alpha2lim)
        Ftil{k} = F0(h) + alpha1lim(i)*F1 + alpha2lim(j)*F2;
        Gtil{k} = G0 + alpha1lim(i)*G1 + alpha2lim(j)*G2;
        CLtil{k}= Ftil{k} - Gtil{k}*K;
        k = k+1;
    end
end
n=3;
cvx_begin sdp
variable M(n,n) symmetric
for i = 1:numel(CLtil)
    [-M, M*CLtil{i}'; CLtil{i}*M, -M] <= -eye(2*n);
end
cvx_end
if strcmp(cvx_status,'Solved')
    P = inv(M);
    
    % Check Positive Definiteness of P and Q
    isposdef = all(eig(P) > 0)
    
    % Check for negative definiteness of  CL'*P*CL - P
    isnegdef = [all((eig(CLtil{1}'*P*CLtil{1} - P)) < 0), ...
        all((eig(CLtil{2}'*P*CLtil{2} - P)) < 0), ...
        all((eig(CLtil{3}'*P*CLtil{3} - P)) < 0), ...
        all((eig(CLtil{4}'*P*CLtil{4} - P)) < 0)]
end

% Gonna get a very similar, though slightly smaller graph to the one in A1.

%%
% K = [22.4 7 0];
hvec = 0.01:0.01:0.47;

% These values are roughly taken from the plot in Assingment 1, and impose
% limits on what intervals to solve for, so we're not just randomly solving
% for every possible tauM and taum combination.
tmin = [zeros(1,39),0.005*ones(1,5),0.01*ones(1,3)];
tmax = [hvec(1:11),0.1*ones(1,4),0.08*ones(1,5),0.06*ones(1,7),0.04*ones(1,14),0.02*ones(1,6)];

% Store the maximum and minimum values of tau for which the polytopic
% overapproximation has a feasible solution.
% tauMvec = zeros(1,numel(hvec));
% taumvec = zeros(1,numel(hvec));
for l = 1:numel(hvec)
    h = hvec(l);
    taum = tmin(l); tauM = tmax(l);
    while 1

        alphasq= [alpha1f(h,tauM), alpha2f(h,tauM);...
                  alpha1f(h,taum), alpha2f(h,tauM);...
                  alpha1f(h,taum), alpha2f(h,taum);
                  alpha1f(h,tauM), alpha2f(h,taum)];
        clear CL
        Ftil{length(alphasq(:,1))} = [];
        Gtil{length(alphasq(:,1))} = [];
        CLtil{length(alphasq(:,1))} = [];
        for k=1:length(alphasq(:,1))
            a1 = alphasq(k,1); a2 = alphasq(k,2);
            Ftil{k} = F0(h) + a1*F1 + a2*F2;
            Gtil{k} = G0 + a1*G1 + a2*G2;
            CLtil{k}= Ftil{k} - Gtil{k}*K;
        end
        n=3;
        cvx_begin sdp quiet
        variable M(n,n) symmetric
        for i = 1:numel(CLtil)
            [-M, M*CLtil{i}'; CLtil{i}*M, -M] <= -eye(2*n);
        end
        cvx_end
        % Check if solved
        if strcmp(cvx_status,'Solved')
            P = inv(M);     
            % Check Positive Definiteness of P and Q
            isposdef = all(eig(P) > 0); 
            % Check for negative definiteness of  CL'*P*CL - P
            isnegdef = [all((eig(CLtil{1}'*P*CLtil{1} - P)) < 0), ...
                all((eig(CLtil{2}'*P*CLtil{2} - P)) < 0), ...
                all((eig(CLtil{3}'*P*CLtil{3} - P)) < 0), ...
                all((eig(CLtil{4}'*P*CLtil{4} - P)) < 0)];
            % Check if solution is just a numerical issue
            if isposdef && all(isnegdef)
                % Store if everything seeems good
                [h tauM]
                tauMvec(l) = tauM;
                taumvec(l) = taum;
                break
            end
        % If not solved
        else
            % Try lowering the maximum tau value
            tauM = tauM - 0.001
            % If the limits are the same, reset tauM and increase taum (to
            % handle h > 0.35).
            if tauM <= taum
                taum = taum+0.001;
                tauM = tmax(l);
            end
            % If taum is the same as the maximum tabled tau, then the thing
            % will have no solution.
            if taum >= tmax(l)
                tauMvec(l) = 0;
                taumvec(l) = 0;
                break
            end
        end
    end
end
%%
load tauMAXvec
load tauMINvec
load goodpts
hvec = 0.01:0.01:0.48;
tauMAXvec = [tauMAXvec 0 0 0];
tauMINvec = [tauMINvec 0 0 0];
% scatter(rows*0.01, (cols-1)*0.001, [], [150 150 150]/256, 'x','LineWidth',2); hold on;
[rows, cols] = find(goodpts == 1);
% Making plots fancy is a fucking art.
hold on;
% plot([hvec(1);hvec(1)], [tauMINvec(1);tauMAXvec(1)],'Color', '#999999', 'Linewidth', 1.5);  
plot([hvec(1:end-3);hvec(1:end-3)], [tauMINvec(1:end-3);tauMAXvec(1:end-3)],'r' , 'Linewidth', 1.5);  
for i = 1:max(rows)-2
    idx = find(rows == i);
    colvec = cols(idx);
    find(colvec == max(colvec));
    if (max(colvec)-1)*0.001 >= tauMAXvec(i) + 10^(-4)
        scatter(i*0.01, (max(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
        plot([i*0.01;i*0.01], [tauMAXvec(i);(max(colvec)-1)*0.001], 'Color','#999999',  'Linewidth', 1.5);
    end
end
for i = max(rows)-1:max(rows)
    idx = find(rows == i);
    colvec = cols(idx);
    find(colvec == max(colvec));
    scatter(i*0.01, (max(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
    scatter(i*0.01, (min(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
    plot([i*0.01;i*0.01], [(min(colvec)-1)*0.001;(max(colvec)-1)*0.001], 'Color','#999999',  'Linewidth', 1.5);
end
scatter(hvec(1:end-3), tauMAXvec(1:end-3), [], 'rx','LineWidth',2); 
scatter(hvec(1:end-3), tauMINvec(1:end-3), [], 'rx','LineWidth',2);
hold off;
grid on; grid minor;
set(gcf,'Position',[10 50 1210 450]);
set(gca, 'Fontsize', 15)
title('Stable combinations of sampling time {\it\fontname{Cambria Math}h} and delay \tau.', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
legend('Stable via Polytopic Over-Approximation', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
%% Question 4.3
K = [22.4 7 0];
hvec = 0.01:0.01:0.47;
mvec = [linspace(0.02,0.105,10), linspace(0.11,0.01,37)];

% These values are roughly taken from the plot in Assingment 1, and impose
% limits on what intervals to solve for, so we're not just randomly solving
% for every possible tauM and taum combination.
tmin = [zeros(1,39),0.005*ones(1,5),0.01*ones(1,3)];
tmax = [hvec(1:11),0.1*ones(1,4),0.08*ones(1,5),0.06*ones(1,7),0.04*ones(1,14),0.02*ones(1,6)];

% Store the maximum and minimum values of tau for which the polytopic
% overapproximation has a feasible solution.
tauMAXvec2 = zeros(1,numel(hvec));
tauMINvec2 = zeros(1,numel(hvec));
for l = 1:length(hvec)
    h = hvec(l);
    taum = tmin(l); tauM = tmax(l);
    while 1
        m = mvec(l);
        f = 10;
        alphahex = [alpha1f(h,taum), alpha2f(h,taum);
            (alpha1f(h,taum)*(1-m/f)+alpha1f(h,tauM)*m/f), alpha2f(h,taum);...
            alpha1f(h,tauM), (alpha2f(h,taum)*m/f+alpha2f(h,tauM)*(1-m/f));...
            alpha1f(h,tauM), alpha2f(h,tauM);...
            (alpha1f(h,taum)*m+alpha1f(h,tauM)*(1-m)), alpha2f(h,tauM);...
            alpha1f(h,taum), (alpha2f(h,taum)*(1-m)+alpha2f(h,tauM)*m)];
        
        clear CLtil Ftil Gtil
        Ftil{size(alphahex,1)} = [];
        Gtil{size(alphahex,1)} = [];
        CLtil{size(alphahex,1)} = [];
        for k=1:size(alphahex,1)
            a1 = alphahex(k,1); a2 = alphahex(k,2);
            Ftil{k} = F0(h) + a1*F1 + a2*F2;
            Gtil{k} = G0 + a1*G1 + a2*G2;
            CLtil{k}= Ftil{k} - Gtil{k}*K;
        end
        n=3;
        cvx_begin sdp quiet
        variable M(n,n) symmetric
        % variable P(n,n) symmetric
        for i = 1:numel(CLtil)
            [-M, M*CLtil{i}'; CLtil{i}*M, -M] <= -eye(2*n);
            M >= eye(n);
            % CLtil{i}'*P*CLtil{i} - P <= -eye(n);
            % P >= eye(n);
        end
        cvx_end
        % Check if solved
        if strcmp(cvx_status,'Solved')
            P = inv(M);     
            % Check Positive Definiteness of P and Q
            isposdef = all(eig(P) > 0); 
            % Check for negative definiteness of  CL'*P*CL - P
            isnegdef = [all((eig(CLtil{1}'*P*CLtil{1} - P)) < 0), ...
                all((eig(CLtil{2}'*P*CLtil{2} - P)) < 0), ...
                all((eig(CLtil{3}'*P*CLtil{3} - P)) < 0), ...
                all((eig(CLtil{4}'*P*CLtil{4} - P)) < 0), ...
                all((eig(CLtil{5}'*P*CLtil{5} - P)) < 0), ...
                all((eig(CLtil{6}'*P*CLtil{6} - P)) < 0)];
            % Check if solution is just a numerical issue
            if isposdef && all(isnegdef)
                % Store if everything seeems good
                [h tauM]
                tauMAXvec2(l) = tauM;
                tauMINvec2(l) = taum;
                break
            end
        % If not solved
        else
            % Try lowering the maximum tau value
            tauM = tauM - 0.001;
            % If the limits are the same, reset tauM and increase taum (to
            % handle h > 0.35).
            if tauM <= taum
                taum = taum+0.001;
                tauM = tmax(l);
            end
            % If taum is the same as the maximum tabled tau, then the thing
            % will have no solution.
            if taum >= tmax(l)
                tauMAXvec2(l) = 0;
                tauMINvec2(l) = 0;
                break
            end
        end
    end
end

save tauMAXvec2
save tauMINvec2
%%
load tauMAXvec2
load tauMINvec2

hvec = 0.01:0.01:0.48;
tauMAXvec2 = [tauMAXvec2 0];
tauMINvec2 = [tauMINvec2 0];
% scatter(rows*0.01, (cols-1)*0.001, [], [150 150 150]/256, 'x','LineWidth',2); hold on;
[rows, cols] = find(goodpts == 1);
% Making plots fancy is a fucking art.
hold on;
% plot([hvec(1);hvec(1)], [tauMINvec(1);tauMAXvec(1)],'Color', '#999999', 'Linewidth', 1.5);  
plot([hvec(1:end-3);hvec(1:end-3)], [tauMINvec2(1:end-3);tauMAXvec2(1:end-3)],'r' , 'Linewidth', 1.5);  
% for i = 1:max(rows)-2
%     idx = find(rows == i);
%     colvec = cols(idx);
%     find(colvec == max(colvec));
%     if (max(colvec)-1)*0.001 >= tauMAXvec(i) + 10^(-4)
%         scatter(i*0.01, (max(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
%         plot([i*0.01;i*0.01], [tauMAXvec(i);(max(colvec)-1)*0.001], 'Color','#999999',  'Linewidth', 1.5);
%     end
% end
% for i = max(rows)-1:max(rows)
%     idx = find(rows == i);
%     colvec = cols(idx);
%     find(colvec == max(colvec));
%     scatter(i*0.01, (max(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
%     scatter(i*0.01, (min(colvec)-1)*0.001, [],  [150 150 150]/256, 'x','LineWidth',2);
%     plot([i*0.01;i*0.01], [(min(colvec)-1)*0.001;(max(colvec)-1)*0.001], 'Color','#999999',  'Linewidth', 1.5);
% end
scatter(hvec(1:end-3), tauMAXvec2(1:end-3), [], 'rx','LineWidth',2); 
scatter(hvec(1:end-3), tauMINvec2(1:end-3), [], 'rx','LineWidth',2);
hold off;
grid on; grid minor;
set(gcf,'Position',[10 50 1210 450]);
set(gca, 'Fontsize', 15)
title('Stable combinations of sampling time {\it\fontname{Cambria Math}h} and delay \tau.', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
legend('Stable via Refined Polytopic Over-Approximation', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
%% Question 5.1
SN = [5 3 2 1 2 6 3];
a = SN(1); b = SN(3); c = SN(end);
A = [a b+0.5; 0 -c];

B = [0; 1];

P = [-2,-3];
K = place(A,B,P); % Find controller with the Place Command

% sigma = 0.5
n=2;
cvx_begin sdp quiet
variable P(n,n) symmetric semidefinite
variable Q(n,n) symmetric semidefinite
% (A-B*K)'*P + P*(A-B*K) >= -Q;
% (A-B*K)'*P + P*(A-B*K) <= -Q;
(A-B*K)'*P + P*(A-B*K) == -Q; 
P >= eye(n);
cvx_end
clc

% Q = r*eye(n); eig(Q)
% Q = 100*eye(n);
isgood = [all(eig(P))>0, all(eig((A-B*K)'*P + P*(A-B*K))<0)]
%% Question 5.2
tstart = 0;
tfinal = 3;
x0 = (rand(2,1)*10)-5; e0 = [0; 0]; 

% x0  = [-1;2.25];
x0  = [2.5;-2.5];
% Condition 1
% xe0 = [x0; x0];

% Condition 2
xe0 = [x0; e0];

sigma = 0.999;
refine = 4;
event = @(t, xe) condition(t, xe, A, B, K, P, Q, sigma);
options = odeset('Events',event,'OutputSel',1,'Refine',refine);

tout = tstart;
xeout = xe0.';
teout = [];
xeeout = [];
ieout = [];

fun = @(t, xe) continuousSystem(t,xe, A, B, K);
while 1
    % Solve until the first terminal event.
    [t,xe,te,xee,ie] = ode45(fun,[tstart tfinal],xe0,options);

    % Accumulate output.  This could be passed out as output arguments.
    nt = length(t);
    tout = [tout; t(2:nt)];
    xeout = [xeout; xe(2:nt,:)];
    teout = [teout; te];          % Events at tstart are never reported.
    xeeout = [xeeout; xee];
    ieout = [ieout; ie];

    if t(end) == tfinal
        break
    end
    % Condition 1
    % xe0 = [xee(1); xee(2); xee(1); xee(2)];

    % Condition 2
    xe0 = [xee(1); xee(2); 0; 0];

    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options = odeset(options,'Events',event,'InitialStep',t(nt)-t(nt-refine),...
        'MaxStep',t(nt)-t(1));


    tstart = t(nt);
end
%%

plot(tout,xeout(:,1), 'Color','#EE0000', 'Linewidth', 1.5); hold on;
plot(tout,xeout(:,2), 'Color','#000000',  'Linewidth', 1.5);
% xline(tLyap, 'k--')
scatter(teout,xeeout(:,1), [], [238,0,0]/256, 'x', 'Linewidth', 2);
scatter(teout,xeeout(:,2),[], [0,0,0]/256, 'x', 'Linewidth', 2); hold off;
grid on; grid minor;
set(gcf,'Position',[10 50 1210 450]);
set(gca, 'Fontsize', 15)
title(['x_0 = [',num2str(x0(1)),',',num2str(x0(2)),...
    '],   \sigma = ',num2str(sigma) ,',   Events: ',num2str(numel(xeeout(:,1)))], ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('Amplitude', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('Time [s]', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
legend('x1', 'x2', 'V(T ) = 0.1V(0)', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');

%% Question 5.3
% P = inv(M);
sims = 200;
X0 = zeros(2,sims);
hvec = zeros(1,sims);
hvec2 = zeros(1,sims);
for j = 1:sims
    
    tstart = 0;
    tfinal = 2;
    x0 = (rand(2,1)*10)-5; 
    % x0 = [-0.25; 1.5];
    e0 = [0; 0]; xe0 = [x0; e0];
    X0(:,j) = x0;
    V0 = x0'*P*x0;
    
    sigma = .5;
    refine = 4;
    event = @(t, xe) condition_Lyap(t, xe,A, B, K, P, Q, sigma, V0);
    options = odeset('Events',event,'OutputSel',1,'Refine',refine);
    
    tout = tstart;
    xeout = xe0.';
    teout = [];
    xeeout = [];
    ieout = [];
    
    fun = @(t, xe) continuousSystem(t,xe, A, B, K);
    flag = 0; % First Time V(T) <= 0.1*V(0)
    while 1
        % Solve until the first terminal event.
        [t,xe,te,xee,ie] = ode45(fun,[tstart,tfinal],xe0,options);
        
        % Accumulate output.  This could be passed out as output arguments.
        nt = length(t);
        tout = [tout; t(2:nt)];
        xeout = [xeout; xe(2:nt,:)];
        teout = [teout; te];          % Events at tstart are never reported.
        xeeout = [xeeout; xee];
        ieout = [ieout; ie];
        
        if t(end) >= tfinal-10^(-2)
            break
        end
        xe0 = [xee(1); xee(2); 0; 0];
        
        if xee(1:2)*P*xee(1:2)' <= 0.1*V0 && flag == 0
            flag = 1;
            xe0 = xee;
            tidx = numel(teout);
            tLyap = te;
            % No longer need to check for Lyapunov decrease
            event = @(t, xe) condition(t, xe,A, B, K, P, Q, sigma);
        end
        % A good guess of a valid first timestep is the length of the last valid
        % timestep, so use it for faster computation.  'refine' is 4 by default.
        options = odeset(options,'Events',event,'InitialStep',t(nt)-t(nt-refine),...
            'MaxStep',t(nt)-t(1));
        
        tstart = t(nt);       
    end
    % Delete event related to Lyapunov Decrease
    hvec(j) = teout(tidx)/(tidx-1);
    if isinf(hvec(j))
        hvec(j) = 0;      
    end
    hvec2(j) = teout(tidx+1)/(tidx);

end
havg = mean(hvec(hvec>0))
havg2 = mean(hvec2(hvec>0))
%%

SYS= ss(A,B,eye(2),[0;0]);
X0 = x0;
T= zeros((100-1)*ceil(tfinal/havg),1); 
X = zeros((100-1)*ceil(tfinal/havg),2);
for i = 1:ceil(tfinal/havg)
t =linspace(0,havg,100); u = -K*X0*ones(1,length(t));
[~,t,x] = lsim(SYS,u,t,X0);
T(99*(i-1)+1:(100-1)*i) = havg*(i-1)+t(1:end-1);
X(99*(i-1)+1:(100-1)*i,:) =  x(1:end-1,:);
X0 = x(end,:)';
end
X(T>tfinal,:) = []; T(T>tfinal) = [];
close 

subplot(2,1,1)
plot(T,X(:,1),'Color', '#EE0000'); hold on;
plot(tout,xeout(:,1),'k');
plot(teout,xeeout(:,1),'ko');
legend('x1 - Discrete', 'x1 - Triggered', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southwest');
ylabel('Amplitude', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on
grid minor
subplot(2,1,2)
plot(T,X(:,2),'Color', '#EE0000'); hold on;
plot(tout,xeout(:,2),'k');
plot(teout,xeeout(:,2),'ko'); hold off;
ylabel('Amplitude', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('Time [s]', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
legend('x2 - Discrete', 'x2 - Triggered', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southwest');
sgtitle(['x_0 = [',num2str(x0(1)),',',num2str(x0(2)),...
    '],   \sigma = ',num2str(sigma) ,',   Events: ',num2str(numel(xeeout(:,1))), ...
    ',   h^{avg}: ', num2str(havg),'s'],  'Fontsize', 20, 'Fontweight', 'bold');
grid on
grid minor
%%

[t,xe,te,xee,ie] = ode45(fun,[tstart,tfinal],xe0,options);

% --------------------------------------------------------------------------

function dxedt = continuousSystem(t,xe, A, B, K)
% Condition 1
% dxedt = [A , -B*K; zeros(2,2), zeros(2,2)]*xe;

% Condition 2
dxedt = [A - B*K, -B*K; -(A-B*K), B*K]*xe;

end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = condition_Lyap(t, xe, A, B, K, P, Q, sigma, V0)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
Phi = [(1-sigma).*Q, P*B*K; K'*B'*P, zeros(2,2)];
value = [xe'*Phi*xe; xe(1:2)'*P*xe(1:2)-0.1*V0]; % detect height = 0

% Condition 3
% sigtil = (1-sigma)*min(eig(Q))/(2*norm(P*B*K,'fro'));
% value1 = sigtil*norm(xe(1:2),'fro') - norm(xe(3:4),'fro');
% value = [value1; xe(1:2)'*P*xe(1:2)-0.1*V0]; % detect height = 0
isterminal = [1;1];   % stop the integration
direction = [-1;-1];   % negative direction
end

function [value,isterminal,direction] = condition(t, xe, A, B, K, P, Q, sigma)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

% Condition 1
% Phi = -[A'*P + P*A + sigma*Q, -P*B*K; -(B*K).'*P, zeros(2,2)];
% value = xe'*Phi*xe; % detect height = 0

% Condition 2
Phi = [(1-sigma)*Q, P*B*K; (B*K).'*P, zeros(2,2)];
value = xe'*Phi*xe; % detect height = 0

% Condition 3
% sigtil = (1-sigma)*min(eig(Q))/(2*norm(P*B*K,'fro'));
% value = sigtil*norm(xe(1:2),'fro') - norm(xe(3:4),'fro');
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
end