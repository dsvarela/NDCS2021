%%
clearvars
close all
clc

%% Problem Formulation for Aircraft 1

LTI1.A=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
LTI1.B=[2 0;0 2;3 0;0 3];
LTI1.x0=[-10;10;-1;1];

LTI2.A=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
LTI2.B=[3 0; 0 3; 7 0; 0 7];
LTI2.x0=[10;10;1;1];

LTI3.A=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
LTI3.B=[1 0; 0 1; 1.1 0; 0 1.1];
LTI3.x0=[10;-10;1;-1];

LTI4.A=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
LTI4.B=[6 0;0 6;20 0; 0 20];
LTI4.x0=[-10;-10;-1;-1];

Tfinal=5;
umax=100;
N = 4; 

% Definition of system dimension
dim.nx=4;     %state dimension
dim.nu=2;     %input dimension
dim.N=Tfinal; %horizon

%Definition of quadratic cost function
weight.Q=eye(dim.nx);   %weight on output
weight.R=eye(dim.nu);   %weight on input

% Generation of prediction model 1
predmod1=predmodgenJ(LTI1,dim);            
[H1,h1]=costgenJ(predmod1,dim,LTI1);

% Generation of prediction model 2
predmod2=predmodgenJ(LTI2,dim);            
[H2,h2]=costgenJ(predmod2,dim,LTI2);

% Generation of prediction model 3
predmod3=predmodgenJ(LTI3,dim);            
[H3,h3]=costgenJ(predmod3,dim,LTI3);

% Generation of prediction model 4
predmod4=predmodgenJ(LTI4,dim);            
[H4,h4]=costgenJ(predmod4,dim,LTI4);

%% Constraints

%Constraints of model 1 
b_eq_1 = predmod1.T(dim.nx*(dim.N-1)+1:end,:)*LTI1.x0;
A_eq_1 = predmod1.S(dim.nx*(dim.N-1)+1:end,:);
b_eq_2 = predmod2.T(dim.nx*(dim.N-1)+1:end,:)*LTI2.x0;
A_eq_2 = predmod2.S(dim.nx*(dim.N-1)+1:end,:);
b_eq_3 = predmod3.T(dim.nx*(dim.N-1)+1:end,:)*LTI3.x0;
A_eq_3 = predmod3.S(dim.nx*(dim.N-1)+1:end,:);
b_eq_4 = predmod4.T(dim.nx*(dim.N-1)+1:end,:)*LTI4.x0;
A_eq_4 = predmod4.S(dim.nx*(dim.N-1)+1:end,:);
b_ineq = ones(2*dim.N*dim.nu,1)*umax/Tfinal;
A_ineq = blkdiag([1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1]);

%% Dual problem 
H1_d = blkdiag(H1,zeros(4,4));
H2_d = blkdiag(H2,zeros(4,4));
H3_d = blkdiag(H3,zeros(4,4));
H4_d = blkdiag(H4,zeros(4,4));
A_eq_1_d = [A_eq_1, -eye(4)];
A_eq_2_d = [A_eq_2, -eye(4)];
A_eq_3_d = [A_eq_3, -eye(4)];
A_eq_4_d = [A_eq_4, -eye(4)];
A_ineq_d =[A_ineq, zeros(20,4)];
%% Individual optimization problems
%initialize
mu12 = zeros(4,1);
mu23 = zeros(4,1);
mu34 = zeros(4,1);
mu41 = zeros(4,1);
a = 1;

for i = 1:250
%Problem 1
    h1_d = [h1, mu12'-mu41'];
    cvx_begin quiet
        variable uopt1(dim.nu*dim.N+4) 
        minimize (uopt1'*H1_d*uopt1 + h1_d*uopt1)
        subject to
        A_ineq_d*uopt1 <= b_ineq;
        A_eq_1_d*uopt1 == - b_eq_1;
    cvx_end
    %uopt1 = quadprog(2*H1_d,h1_d,A_ineq_d,b_ineq,A_eq_1_d,-b_eq_1);
%Problem 2
    h2_d = [h2, mu23'-mu12'];
    cvx_begin quiet
        variable uopt2(dim.nu*dim.N+4) 
        minimize (uopt2'*H2_d*uopt2 + h2_d*uopt2)
        subject to
        A_ineq_d*uopt2 <= b_ineq;
        A_eq_2_d*uopt2 == - b_eq_2;
    cvx_end
    %uopt2 = quadprog(2*H2_d,h2_d,A_ineq_d,b_ineq,A_eq_2_d,-b_eq_2);
%Problem 3
    h3_d = [h3, mu34'-mu23'];
    cvx_begin quiet
        variable uopt3(dim.nu*dim.N+4) 
        minimize (uopt3'*H3_d*uopt3 + h3_d*uopt3)
        subject to
        A_ineq_d*uopt3 <= b_ineq;
        A_eq_3_d*uopt3 == - b_eq_3;
    cvx_end
    %uopt3 = quadprog(2*H3_d,h3_d,A_ineq_d,b_ineq,A_eq_3_d,-b_eq_3);
%Problem 4
    h4_d = [h4, mu41'-mu34'];
    cvx_begin quiet
        variable uopt4(dim.nu*dim.N+4) 
        % minimize (uopt4'*H4_d*uopt4 + h4_d*uopt4)
        minimize (uopt4'*H4_d*uopt4 + h4_d*uopt4)
        subject to
        A_ineq_d*uopt4 <= b_ineq;
        A_eq_4_d*uopt4 == - b_eq_4;
    cvx_end
    %uopt4 = quadprog(2*H4_d,h4_d,A_ineq_d,b_ineq,A_eq_4_d,-b_eq_4);
    %mu update
    mu12 = mu12 + a*(uopt1(end-3:end)-uopt2(end-3:end));
    mu23 = mu23 + a*(uopt2(end-3:end)-uopt3(end-3:end));
    mu34 = mu34 + a*(uopt3(end-3:end)-uopt4(end-3:end));
    mu41 = mu41 + a*(uopt4(end-3:end)-uopt1(end-3:end));
    mu_save(:,i) = [mu12;mu23;mu34;mu41];
end

%% Accelerated gradient 

L = norm(A_eq_1*H1*A_eq_1',2);
%% Centralized solution 

Hc = blkdiag(H1,H2,H3,H4,zeros(4,4));
hc = [h1, h2, h3, h4,zeros(1,4)];
A_eq_c1 = [blkdiag(A_eq_1,A_eq_2,A_eq_3,A_eq_4)];
A_eq_c2 = [-eye(4);-eye(4);-eye(4);-eye(4)];
A_eq_c = [A_eq_c1, A_eq_c2];
b_eq_c = [b_eq_1;b_eq_2;b_eq_3;b_eq_4];
A_ineq_c = [blkdiag(A_ineq,A_ineq,A_ineq,A_ineq),zeros(80,4)];
b_ineq_c = ones(2*4*dim.nu*dim.N,1)*umax/Tfinal;

cvx_begin
    variable uopt(4*dim.nu*dim.N+4) 
    minimize (uopt'*Hc*uopt + hc*uopt)
    subject to
    A_ineq_c*uopt <= b_ineq_c;
    A_eq_c*uopt == -b_eq_c;
cvx_end

%% Aircarft state trajectories for centralized case

u1 = uopt(1:10);
x1 = [LTI1.x0; predmod1.T*LTI1.x0+predmod1.S*u1];
x11 = x1(1:4:end);
x12 = x1(2:4:end);
x13 = x1(3:4:end);
x14 = x1(4:4:end);
u2 = uopt(11:20);
x2 = [LTI2.x0;predmod2.T*LTI2.x0+predmod2.S*u2];
x21 = x2(1:4:end);
x22 = x2(2:4:end);
x23 = x2(3:4:end);
x24 = x2(4:4:end);
u3 = uopt(21:30);
x3 = [LTI3.x0;predmod3.T*LTI3.x0+predmod3.S*u3];
x31 = x3(1:4:end);
x32 = x3(2:4:end);
x33 = x3(3:4:end);
x34 = x3(4:4:end);
u4 = uopt(31:40);
x4 = [LTI4.x0; predmod4.T*LTI4.x0+predmod4.S*u4];
x41 = x4(1:4:end);
x42 = x4(2:4:end);
x43 = x4(3:4:end);
x44 = x4(4:4:end);
xf = uopt(end-3:end);
T = [0,1,2,3,4,5];

% Plot the trajectories
figure(1) 
hold on 
plot(T,x11,'r');
plot(T,x21,'r:');
plot(T,x31,'r--');
plot(T,x41,'r-.');
plot(T,x12,'b');
plot(T,x22,'b:');
plot(T,x32,'b--');
plot(T,x42,'b-.');
plot(T,x13,'g');
plot(T,x23,'g:');
plot(T,x33,'g--');
plot(T,x43,'g-.');
plot(T,x14,'c');
plot(T,x24,'c:');
plot(T,x34,'c--');
plot(T,x44,'c-.');

%% ADMM

%initialize
xf0 = [-5;-3;-7;-4]; %optimal solution from combined -5.9566   -3.7671   -7.3848   -4.6138
y0 = [0;0;0;0];
rho = 10;
iter = 200;

xf = zeros(N,iter);
xf(:,1) = xf0;
y1 = zeros(dim.nx, iter);
y1(:,1) = y0;
y2 = zeros(dim.nx, iter);
y2(:,1) = y0;
y3 = zeros(dim.nx, iter);
y3(:,1) = y0;
y4 = zeros(dim.nx, iter);
y4(:,1) = y0;
uopt1_ = zeros(iter, dim.nu*dim.N);
uopt2_ = zeros(iter, dim.nu*dim.N);
uopt3_ = zeros(iter, dim.nu*dim.N);
uopt4_ = zeros(iter, dim.nu*dim.N);
ubar = zeros(N,iter);
ubar(:,1) = y0;
%% ADMM algorithm
for k = 1:iter
    
    %1) Update x's
    %x1(k+1)
    cvx_begin
            variable uopt1(dim.nu*dim.N) 
            minimize (uopt1'*H1*uopt1 + h1*uopt1 + y1(:,k)'*(A_eq_1*uopt1 + b_eq_1 - xf(:,k)) + rho/2*sum_square(A_eq_1*uopt1 + b_eq_1 - xf(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
            subject to
            for i = 0:dim.N-1
                norm(uopt1(i*2+1:i*2+2),2) <= umax/Tfinal;
            end
    cvx_end
    uopt1_(k,:) = uopt1;
    %x2(k+1)
    cvx_begin
            variable uopt2(dim.nu*dim.N) 
            minimize (uopt2'*H2*uopt2 + h2*uopt2 + y2(:,k)'*(A_eq_2*uopt2 + b_eq_2 - xf(:,k)) + rho/2*sum_square(A_eq_2*uopt2 + b_eq_2 - xf(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
            subject to
            for i = 0:dim.N-1
                norm(uopt2(i*2+1:i*2+2),2) <= umax/Tfinal;
            end
    cvx_end
    uopt2_(k,:) = uopt2;
    %x3(k+1)
    cvx_begin
            variable uopt3(dim.nu*dim.N) 
            minimize (uopt3'*H3*uopt3 + h3*uopt3 + y3(:,k)'*(A_eq_3*uopt3 + b_eq_3 - xf(:,k)) + rho/2*sum_square(A_eq_3*uopt3 + b_eq_3 - xf(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
            subject to
            for i = 0:dim.N-1
                norm(uopt3(i*2+1:i*2+2),2) <= umax/Tfinal;
            end
    cvx_end
    uopt3_(k,:) = uopt3;
    %x4(k+1)
    cvx_begin
            variable uopt4(dim.nu*dim.N) 
            minimize (uopt4'*H4*uopt4 + h4*uopt4 + y4(:,k)'*(A_eq_4*uopt4 + b_eq_4 - xf(:,k)) + rho/2*sum_square(A_eq_4*uopt4 + b_eq_4 - xf(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
            subject to
            for i = 0:dim.N-1
                norm(uopt4(i*2+1:i*2+2),2) <= umax/Tfinal;
            end
    cvx_end
    uopt4_(k,:) = uopt4;
    
    %2) Update xf
    xf(:,k+1) = 1/N*((A_eq_1*uopt1 + b_eq_1 + 1/rho.*y1(:,k))+(A_eq_2*uopt2 + b_eq_2 + 1/rho.*y2(:,k))+(A_eq_3*uopt3 + b_eq_3 + 1/rho.*y3(:,k))+(A_eq_4*uopt4 + b_eq_4 + 1/rho.*y4(:,k)));
    
    %3) Update y's
    y1(:,k+1) = y1(:,k) + rho*(A_eq_1*uopt1 + b_eq_1 - xf(:,k+1));
    y2(:,k+1) = y2(:,k) + rho*(A_eq_2*uopt2 + b_eq_2 - xf(:,k+1));
    y3(:,k+1) = y3(:,k) + rho*(A_eq_3*uopt3 + b_eq_3 - xf(:,k+1));
    y4(:,k+1) = y4(:,k) + rho*(A_eq_4*uopt4 + b_eq_4 - xf(:,k+1));
end

%% Second version of ADMM algorithm
% for k = 1:iter
%     
%     %1) Update x's
%     %x1(k+1)
%     cvx_begin
%             variable uopt1(dim.nu*dim.N) 
%             minimize (uopt1'*H1*uopt1 + h1*uopt1 + y1(:,k)'*(A_eq_1*uopt1 + b_eq_1 - ubar(:,k)) + rho/2*sum_square(A_eq_1*uopt1 + b_eq_1 - ubar(:,k)))
%             subject to
%             for i = 0:dim.N-1
%                 norm(uopt1(i*2+1:i*2+2),2) <= umax/Tfinal;
%             end
%     cvx_end
%     uopt1_(k,:) = uopt1;
%     %x2(k+1)
%     cvx_begin
%             variable uopt2(dim.nu*dim.N) 
%             minimize (uopt2'*H2*uopt2 + h2*uopt2 + y2(:,k)'*(A_eq_2*uopt2 + b_eq_2 - ubar(:,k)) + rho/2*sum_square(A_eq_2*uopt2 + b_eq_2 - ubar(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
%             subject to
%             for i = 0:dim.N-1
%                 norm(uopt2(i*2+1:i*2+2),2) <= umax/Tfinal;
%             end
%     cvx_end
%     uopt2_(k,:) = uopt2;
%     %x3(k+1)
%     cvx_begin
%             variable uopt3(dim.nu*dim.N) 
%             minimize (uopt3'*H3*uopt3 + h3*uopt3 + y3(:,k)'*(A_eq_3*uopt3 + b_eq_3 - ubar(:,k)) + rho/2*sum_square(A_eq_3*uopt3 + b_eq_3 - ubar(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
%             subject to
%             for i = 0:dim.N-1
%                 norm(uopt3(i*2+1:i*2+2),2) <= umax/Tfinal;
%             end
%     cvx_end
%     uopt3_(k,:) = uopt3;
%     %x4(k+1)
%     cvx_begin
%             variable uopt4(dim.nu*dim.N) 
%             minimize (uopt4'*H4*uopt4 + h4*uopt4 + y4(:,k)'*(A_eq_4*uopt4 + b_eq_4 - ubar(:,k)) + rho/2*sum_square(A_eq_4*uopt4 + b_eq_4 - ubar(:,k))) % - lambda'*(A_eq_1*uopt1 + b_eq_1 - xf))
%             subject to
%             for i = 0:dim.N-1
%                 norm(uopt4(i*2+1:i*2+2),2) <= umax/Tfinal;
%             end
%     cvx_end
%     uopt4_(k,:) = uopt4;
%     
%     %2) Update xbar
%     ubar(:,k+1) = 1/N*(A_eq_1*uopt1 + b_eq_1 + A_eq_2*uopt2 + b_eq_2 + A_eq_3*uopt3 + b_eq_3 + A_eq_4*uopt4 + b_eq_4);
%     
%     %3) Update y's
%     y1(:,k+1) = y1(:,k) + rho*(A_eq_1*uopt1 + b_eq_1 - ubar(:,k+1));
%     y2(:,k+1) = y2(:,k) + rho*(A_eq_2*uopt2 + b_eq_2 - ubar(:,k+1));
%     y3(:,k+1) = y3(:,k) + rho*(A_eq_3*uopt3 + b_eq_3 - ubar(:,k+1));
%     y4(:,k+1) = y4(:,k) + rho*(A_eq_4*uopt4 + b_eq_4 - ubar(:,k+1));
% end
% xf = A_eq_1*uopt1 + b_eq_1

%% plot results 

figure(3)
plot(xf(1,:))
hold on 
plot(xf(2,:))
plot(xf(3,:))
plot(xf(4,:))

%% Error sequence plot 
e_xf = zeros(4,1);
for i = 1:size(xf,2)-1;
    e_xf(:,i) = xf(:,i+1)-xf(:,i);
end

figure(4)
plot(e_xf(1,:))
hold on
plot(e_xf(2,:))
plot(e_xf(3,:))
plot(e_xf(4,:))