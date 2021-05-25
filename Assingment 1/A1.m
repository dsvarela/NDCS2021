clearvars; clc; close all;
%% Student Number
SN = [5 3 2 1 2 6 3];

a = SN(1); b = SN(3); c = SN(end);
syms h s t
%% Matrices
A = [a b+0.5; 0 -c];
B = [0; 1];

%% Question 1.1
P = [-2,-3];
K = place(A,B,P); % Find controller with the Place Command
AK = A-B*K; % Close the loop

%% Question 1.2 (calculate the matrices and close the loop)
% Using sybolic variables to calculate exact expressions

% This whole thing checks out, and the derivation seems correct.
Fx = expm(A*h); % Fx = e^(Ah)
G = int(expm(A*s)*[0;1],s,0,h); % G = int_{0}^{h} e^(As) ds

FxK = Fx-G*K; % Closed loop
eig_h = eig(FxK); % Find expressions for eigenvalues of the closed loop.

%% Question 1.3 (Plotting Things)
close all
figure(1);
set(gcf,'Position',[10 50 1210 450]);
fplot(abs(eig_h(1)),[0,0.4],'-', 'Color','#aaaaaa', 'Linewidth', 1.5);
hold on
fplot(abs(eig_h(2)),[0,0.4],'-', 'Color','#000000', 'Linewidth', 1.5);
fplot(max(abs(eig_h(1)),abs(eig_h(2))),[0,0.4],'-', 'Color','#EE0000', 'Linewidth', 1.5);
legend('|\lambda_1|','|\lambda_2|','max |\lambda_i|', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southwest');
set(gca, 'Fontsize', 15)
title('Evolution of Eigenvalues \lambda with Sampling Time {\it\fontname{Cambria Math}h}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('|\lambda|', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')

grid on
grid minor
%% Question 2.1
% Turn them into functions so it's easier to keep track of what I'm
% actually doing here.

M0  = int(expm(A*s)*[0;1],s,0,h-t); % M0 = int_{0}^{h-tau} e^(As) ds (G1)
M1  = int(expm(A*s)*[0;1],s,h-t,h); % M1 = int_{h-tau}^{h} e^(As) ds (Fu)

% Extended system for small delays
Fe = [[Fx, M1]; zeros(1,3)];
Ge = [M0; 1];
FKe = Fe - Ge*[K 0]; % expressions for eigenvalues of the extended closed loop.

%  Turn eigenvalues expression into a function we can evaluate
eigf_Ke = matlabFunction(eig(FKe)); 

% Loop for values of h and tau, calculate eigf_Ke(h,t) and store the larges
% eigenvalue for each pair (h,tau)
hvec = 0.01:0.01:0.5;
maxeig = zeros(length(hvec), length(0:0.005:hvec(end)));
goodpts =  zeros(length(hvec), length(0:0.005:hvec(end)));
for i = 1:length(hvec)
    hn = hvec(i);
    % Crop at 0.125 cause I know combinations higher than that never
    % yield stable results (I know cause I tried it and it didn't work)
    tvec = 0:0.001:min(hn,0.125+eps);
    for j = 1:length(tvec)
        tn = tvec(j);
        maxeig(i,j) = max(abs(eigf_Ke(hn,tn)));
        if maxeig(i,j) < 1
            goodpts(i,j) = 1;
        end
    end
end

% Find rows and columns for which the eigenvalues are smaller than 1
[rows,cols,~] = find(goodpts == 1);

% Plot Things
figure(2);
set(gcf,'Position',[10 50 1210 450]);

c = [238, 0, 0]/256;
scatter(rows*0.01, cols*0.001-0.001, [],  [238, 0, 0]/256, 'x');
% [rows,cols,~] = find(maxeig > 1);
% scatter(rows*0.01, cols*0.001, [],  [170, 170, 170]/256, 'x');

legend('max |\lambda| < 1', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Stable Combinations of Sampling Time {\it\fontname{Cambria Math}h} and Delay \tau', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on; grid minor;
%% Question 2.3 - Dynamic Controller + LQR

%  Turn symbolic expressions functions we can evaluate
Fxf = matlabFunction(Fx);
Gf = matlabFunction(G);
M0f = matlabFunction(M0);
M1f = matlabFunction(M1);

hn = 0.2; % Fixed h, picked +/- at random from the range of stable values
Q = 0.1*eye(2); R = 10; % Tuned iteratively
% Find the first two elements of K = [k1 k2 k3]. These correspond to the
% values of K found for the continuous time.
Kb = dlqr(Fxf(hn), Gf(hn), Q, R);  

% Tune k3 iteratively.
Kb = [Kb 0.35]; % h = 0.2 (+ 0.0878 s)
% Kb = [21.9410    6.8574    0.3612]; % From LMIs
% Decent values of k3 for other sampling times.
% Kb = [Kb 0.4]; % h = 0.175 (+ 0.1011 s)
% Kb = [Kb 0.28]; % h = 0.25 (+ 0.0703 s)
% Kb = [Kb 0.2]; % h = 0.3 (+ 0.0502 s)
% Kb = [Kb 0.17]; % h = 0.35 (+ 0.0428 s)

% Loop for different values of delay, calculate and store the largest
% eigenvalue for each
tvec = 0:0.0001:hn;
maxeig = zeros(5,length(tvec));
for i = 1:length(tvec)
    tn = tvec(i);

    Fe = [[Fxf(hn), M1f(hn,tn)]; zeros(1,3)];
    Ge = [M0f(hn,tn); 1];
    
    maxeig(1,i) = max(abs(eig(Fe-Ge*[K 0])));
    maxeig(2,i) = max(abs(eig(Fe-Ge*Kb)));
end

% Find delay at which the original system crosses |\lambda| = 1
tidx = find(maxeig(1,:) > 0.999 & maxeig(1,:) < 1.001);
maxK = tvec(round(mean(tidx)));
% Find delay at which the modified system crosses |\lambda| = 1
tidx = find(maxeig(2,:) > 0.999 & maxeig(2,:) < 1.001);
maxKb = tvec(round(mean(tidx)));
% Find the difference to evaluate the performance of the controller.
disp(maxKb - maxK);

% Plot Things
figure(3);
set(gcf,'Position',[10 50 1210 450]);

plot(tvec, maxeig(1,:), 'Color', '#AAAAAA','Linewidth', 1.5); hold on;
plot(tvec, maxeig(2,:), 'Color', '#EE0000','Linewidth', 1.5);
plot(tvec, ones(1,length(tvec)), '--k'); hold off
xlim([0,hn]); grid on; grid minor;

legend('K = [K̅ 0]','K = [K_{LQR} 0.35]', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southeast');
set(gca, 'Fontsize', 15)
title('Evolution of Eigenvalues \lambda with Delay \tau for fixed {\it\fontname{Cambria Math}h} = 0.2{\it\fontname{Cambria Math}s}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('|\lambda|', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')

%% Question 2.3 (Solved via LMIs)
hn = 0.2; tvec = 0:0.005:hn;
tvec(end-9:end) = [];
[Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f);
for i = 1:length(Fext)
   Fext{i}(:,end) = [];
   Fext{i}(end,:) = []; 
   Gext{i}(end) = []; 
end
n = 3;
cvx_begin sdp
variable M(n,n) symmetric
variable Y(1,n)
for j = 1:length(Fext)
    [-M, M*Fext{j}'-Y'*Gext{j}'; Fext{j}*M-Gext{j}*Y, -M] <= -eye(2*n);
end
cvx_end
Klmi = Y/M;

%% Question 3.1

% Extended system for longer delays
% The way this matrix is constructed, it accepts values of t \in [0,h),
% even though the value real value is tau \in [h,2*h).
Fe = [[Fx, M0, M1]; 0 0 0 0; 0 0 1 0];
Ge = [0; 0; 1; 0];
FKe = Fe - Ge*[K 0 0]; % expressions for eigenvalues of the extended closed loop.

%  Turn eigenvalues into a function we can evaluate
eigf_Ke = matlabFunction(eig(FKe)); 

% Loop for values of h and tau, calculate eigf_Ke(h,t) and store the larges
% eigenvalue for each pair (h,tau)
maxeig = zeros(length(hvec), length(0:0.005:hvec(end)));
hvec = 0.01:0.01:0.2;
for i = 1:length(hvec)
    hn = hvec(i);
    tvec = 0:0.001:hn;
    for j = 1:length(tvec)
        tn = tvec(j);
        % Crop at 0.125 cause I know combinations higher than that never
        % yield stable results (I know cause I tried it and it didn't work)
        if tn+hn <= 0.125
            maxeig(i,j) = max(abs(eigf_Ke(hn,tn)));
        end
    end
end
% Find rows and columns for which the eigenvalues are smaller than 1
[rows,cols,~] = find(maxeig < 1 & maxeig > 0);

% Plot Things (on the same figure as Question 2.1).
figure(2); hold on;
% Also adding back the h to t before plotting, to get the real delay
% tau = t+h \in [h,2*h)
scatter(rows*0.01, cols*0.001+rows*0.01, [], [0 0 0], 'x'); hold off;
% [rows,cols,~] = find(maxeig > 1);
% scatter(rows*0.01, cols*0.001+rows*0.01, 'rx');

legend('{\fontname{Cambria Math}max |\lambda| < 1 (\tau \leq {\it h})}',...
    '{\fontname{Cambria Math}max |\lambda| < 1 (\tau >{\it h})}',...
    'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northeast');
set(gca, 'Fontsize', 15)
title('Stable Combinations of Sampling Time {\it\fontname{Cambria Math}h} and Delay \tau', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
%% Question 3.3 % Dynamic Controller + LQR

hn = 0.1; % Again picked +/- at random from the set of stable values.
tn = hn; % The same as the hn

Kb = dlqr(Fxf(hn), Gf(hn), 0.1*eye(2), 10);
Kb = [Kb 0.32 0.25]; % h = 0.1 (+ 0.0894 s (max))
Klmi = [37.4095   11.7059    0.6918    0.4706];
% Loop for different values of delay, calculate and store the largest
% eigenvalue for each
tvec = 0:0.0001:hn;
maxeig = zeros(2,length(tvec));
for i = 1:length(tvec)
    tn = tvec(i);
    
    Fe = [[Fxf(hn), M0f(hn,tn), M1f(hn,tn)]; 0 0 0 0; 0 0 1 0];
    Ge = [0; 0; 1; 0]; 
    
    maxeig(1,i) = max(abs(eig(Fe-Ge*[K 0 0])));
    maxeig(2,i) = max(abs(eig(Fe-Ge*Kb)));
    maxeig(3,i) = max(abs(eig(Fe-Ge*Klmi)));
end

% Find delay at which the original system crosses |\lambda| = 1
tidx = find(maxeig(1,:) > 0.999 & maxeig(1,:) < 1.001);
maxK = tvec(round(mean(tidx)));
% Find delay at which the modified system crosses |\lambda| = 1
tidx = find(maxeig(2,:) > 0.999 & maxeig(2,:) < 1.001);
if isnan(round(mean(tidx)))
    tf = length(tvec);
else
    tf = round(mean(tidx));
end
maxKb = tvec(tf);
disp(maxKb - maxK);

% Plot things
figure(4);
set(gcf,'Position',[10 50 1210 450]);

plot(tvec+hn, maxeig(1,:), 'Color', '#AAAAAA','Linewidth', 1.5); hold on;
plot(tvec+hn, maxeig(2,:), 'Color', '#0000EE','Linewidth', 1.5);
plot(tvec+hn, maxeig(3,:), 'Color', '#EE0000','Linewidth', 1.5);
plot(tvec+hn, ones(1,length(tvec)), 'k--'); hold off
xlim([hn,2*hn]); grid on; grid minor;

legend('K = [K̅ 0]','K = [K_{LQR} K_u]', 'K_{LMI}', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'northwest');
set(gca, 'Fontsize', 15)
title('Evolution of Eigenvalues \lambda with Delay \tau for fixed {\it\fontname{Cambria Math}h} = 0.1{\it\fontname{Cambria Math}s}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('|\lambda|', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('\tau', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')

%% Question 3.3 (Solved with LMIs)
hn = 0.1; tvec = hn:0.001:2*hn;
[Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f);

n = 4;
cvx_begin sdp
variable M(n,n) symmetric
variable Y(1,n)
for j = 1:length(Fext)
    [-M, M*Fext{j}'-Y'*Gext{j}'; Fext{j}*M-Gext{j}*Y, -M] <= -eye(2*n);
end
cvx_end
Klmi = Y/M

%% Question 4.1
% K = [22.4 7]; % Original controller

% Very much picked at random (needs to be able to handle 1.5*h of delay)
hn = 0.077;
% hn = 0.106;
tvec = [0.2*hn; 0.5*hn; hn; 1.5*hn]; % delay range

% Loop for the 4 delays and find the respective augmented matrices F and G,
% plus the closed loop CL = F-G*K. All matrices consider xe = [x u_{k-1} u_{k-2}]

[Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f);
clear CL; CL{4}=[];
for i = 1:4
    CL{i} = Fext{i} - Gext{i}*[K 0 0];
end
n = 4;

% This method (Brute Force) works, but matrices are pretty ill-conditioned.
% cvx_begin sdp
% variable P(n,n) symmetric
% CL{1}'*P*CL{1} - P<= -eye(n);
% CL{2}'*P*CL{2} - P<= -eye(n);
% CL{3}'*P*CL{3} - P<= -eye(n);
% CL{4}'*P*CL{4} - P<= -eye(n);
% cvx_end

% Tried it also in matrix form, and it seems to work better
cvx_begin sdp quiet
variable M(n,n) symmetric
[-M, M*CL{1}'; CL{1}*M, -M] <= -eye(2*n);
[-M, M*CL{2}'; CL{2}*M, -M] <= -eye(2*n);
[-M, M*CL{3}'; CL{3}*M, -M] <= -eye(2*n);
[-M, M*CL{4}'; CL{4}*M, -M] <= -eye(2*n);
cvx_end
P = inv(M);

% Check Positive Definiteness of P and Q
isposdef = all(eig(P) > 0)

% Check for negative definiteness of  CL'*P*CL - P
isnegdef = [all((eig(CL{1}'*P*CL{1} - P)) < 0), ...
            all((eig(CL{2}'*P*CL{2} - P)) < 0), ...
            all((eig(CL{3}'*P*CL{3} - P)) < 0), ...
            all((eig(CL{4}'*P*CL{4} - P)) < 0)]

%% Question 4.2 (the LMI way to do it)
% Working, but results are not that good...

n = 4;
hvec = 0.13:0.001:0.135;
for i = 1:length(hvec)
    
    hn = hvec(i);
    tvec = [0.2*hn; 0.5*hn; hn; 1.5*hn];
    [Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f);
    
    cvx_begin sdp
    variable M(n,n) symmetric
    variable Y(1,n)
    % variable Q(2*n,2*n)
    for j = 1:length(Fext)
        [-M, M*Fext{j}'-Y'*Gext{j}'; Fext{j}*M-Gext{j}*Y, -M] <= -eye(2*n,2*n);
    end
    cvx_end
    if ~strcmp(cvx_status,'Infeasible')
        hb = hn;
        Mb = M;
        Yb = Y;
    end
end
hb
P = inv(Mb);
isposdef = all(eig(P) > 0)
Klmi = Yb/Mb;
%% 
hvec = 0.001:0.001:0.35;
maxeig = zeros(2,length(hvec));
for i = 1:length(hvec)
    tvec = [0.2*hvec(i); 0.5*hvec(i); hvec(i)];
    [Fext, Gext] = getmats(hvec(i),tvec,Fxf,M1f,M0f);
    for j = 1:length(tvec)
        CL{j} = Fext{j} - Gext{j}*[Klmi];
        maxeig(j,i) = max(abs(eig(CL{j})));
    end
%         CLi = Fext{end} - Gext{3}*[K 0 0];
%         maxeig(1,i) = max(abs(eig(CLi)));
%         CLi = Fext{end} - Gext{3}*[Klmi];
%         maxeig(2,i) = max(abs(eig(CLi)));
end

tidx = find(maxeig(1,:) < 1.001, 1, 'last');
maxK = hvec(tidx);

tidx = find(maxeig(2,:) < 1.001, 1, 'last');
maxKb = hvec(tidx);
disp(maxKb - maxK);

% plot(hvec, maxeig(1,:)); hold on;
% plot(hvec, maxeig(2,:));
% plot(hvec, maxeig(3,:));
% % plot(hvec, maxeig(4,:));
% % plot(hvec, ones(1,length(hvec))); hold off
% xlim([0,hvec(end)]); grid on; grid minor;

%% Question 4.4

n = 4;
hvec = 0.19:0.001:0.195;
for i = 1:length(hvec)
    
    hn = hvec(i);
    tvec = [0.2*hn; 0.5*hn; hn];
    [Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f);
    
    cvx_begin sdp quiet
    variable M(n,n) symmetric
    variable Y(1,n)

    for j = 1:length(Fext)
        [-M, M*Fext{j}'-Y'*Gext{j}'; Fext{j}*M-Gext{j}*Y, -M] <= -eye(2*n);
    end
    cvx_end
    if ~strcmp(cvx_status,'Infeasible')
        hb = hn;
        Mb = M;
        Yb = Y;
    end
end

P = inv(Mb);
isposdef = all(eig(P) > 0)
Klmi = Yb/Mb;

%% Question 4.3
hvec = 0.001:0.001:0.3
maxeig = zeros(2,length(hvec));
for i = 1:length(hvec)
h=hvec(i);
tvec = [0.2*h; h; 0.5*h];
[Fext, Gext] = getmats(h,tvec,Fxf,M1f,M0f);
for j = 1:length(tvec)
    CL{j} = Fext{j} - Gext{j}*Klmi;
end
CLper = CL{3}*CL{2}*CL{1};
maxeig(1,i) = max(abs(eig(CLper)));

for j = 1:length(tvec)
    CL{j} = Fext{j} - Gext{j}*[22.4, 7 0 0];
end
CLper = CL{3}*CL{2}*CL{1};
maxeig(2,i) = max(abs(eig(CLper)));

end

figure(5);
set(gcf,'Position',[10 50 1210 450]);
plot(hvec,maxeig(1,:), '-', 'Color','#EE0000', 'Linewidth', 1.5); hold on
plot(hvec,maxeig(2,:), '-', 'Color','#AAAAAA', 'Linewidth', 1.5); hold on
ylim([0,1.5])
plot(hvec,ones(length(hvec)), '--', 'Color','#000000', 'Linewidth', 0.1); hold off
legend('max |\lambda_i| (LMI)', 'max |\lambda_i|  ([K̅ 0])', 'Fontweight', 'bold', ...
    'Fontsize', 20, 'Interpreter', 'Tex', 'Location', 'southwest');
set(gca, 'Fontsize', 15)
title('Evolution of Eigenvalues \lambda with Sampling Time {\it\fontname{Cambria Math}h}', ...
    'Fontsize', 20, 'Fontweight', 'bold')
ylabel('|\lambda|', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
xlabel('{\it\fontname{Cambria Math}h}', 'Interpreter', 'Tex',...
    'Fontsize', 20, 'Fontweight', 'bold')
grid on
grid minor
% % Tried it with "simplified LMI", aka just a single one, and it also works
% cvx_begin sdp quiet
% variable M(n,n) symmetric
% [-M, M*CLper'; CLper*M, -M] <= -eye(2*n);
% cvx_end
% P = inv(M);
% % Check Positive Definiteness of P and Q
% isposdef = all(eig(P) > 0)



function [Fext, Gext] = getmats(hn,tvec,Fxf,M1f,M0f)
Fext{length(tvec)}=[]; Gext{length(tvec)}=[]; 
for i = 1:length(tvec)
    tn = tvec(i);
    if tn < hn
        Fext{i} = [Fxf(hn), M1f(hn,tn), [0;0]; 0, 0, 0, 0; 0, 0, 1, 0];
        Gext{i} = [M0f(hn,tn); 1; 0]; 
    else
        tn = tn-hn; % This structure is used for t = [h,2h), so we remove the offset.
        Fext{i} = [Fxf(hn), M0f(hn,tn), M1f(hn,tn); 0, 0, 0, 0; 0, 0, 1, 0];
        Gext{i} = [0; 0; 1; 0];
    end
end
end