clearvars; clc; close all;
%% Student Number
SN = [5 3 2 1 2 6 3];

a = SN(1); b = SN(3); c = SN(end);

%% Matrices
A = [a b+0.5; 0 -c];
B = [0; 1];

%% Question 1.1
P = [-2,-3];
K = place(A,B,[-2,-3]);
AK = A-B*K;

%% Question 1.2
% This whole thing checks out, and the derivation seems correct.
h = 0.1;
Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
Gamma = [(exp(5*h)/16 + 5*exp(-3*h)/48-1/6); (1/3 -exp(-3*h)/3)];
PhiK = Phi-Gamma*K; 
eig_h = eig(PhiK);

%% Question 1.3
hlen = 1000;
hvec = linspace(0.0001, 0.4, hlen);
maxeig = zeros(size(hvec));
K = [22.4 7];
for i = 1:hlen
    h = hvec(i);
    Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
    Gamma = [(exp(5*h)/16 + 5*exp(-3*h)/48-1/6); (1/3 -exp(-3*h)/3)];
    L = eig(Phi - Gamma*K);
    maxeig(i) = max(abs(L));
end

plot(hvec,maxeig);

%% Question 2.1
% Turn them into functions so it's easier to keep track of what I'm
% actually doing here.
%% Question 2.1
% Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
% Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];

tic
K = [22.4 7];
hvec = 0.01:0.01:0.5;
maxeig = zeros(length(hvec), length(0:0.005:hvec(end)));
for i = 1:length(hvec)
    h = hvec(i);
    tvec = 0:0.001:h;
    for j = 1:length(tvec)
        t = tvec(j);
        if t <= 0.125+eps
            Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
            Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
            Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
            Phie = [Phi-Gamma0*K, -Gamma1*K; eye(size(Phi)), zeros(size(Phi))];
            L = eig(Phie); maxeig(i,j) = max(abs(L));
            %           Do this exactly like it's done in the lectures
            %             Phie = [[Phi, Gamma1]; zeros(1,3)];
            %             Gammae = [Gamma0; 1];
            %             L = eig(Phie- Gammae*[K 0]); maxeig(i,j) = max(abs(L));
        end
    end
end
[rows,cols,~] = find(maxeig < 1 & maxeig > 0);
scatter(rows*0.01, cols*0.001, 'gx'); hold on
[rows,cols,~] = find(maxeig > 1);
scatter(rows*0.01, cols*0.001, 'rx');
toc
%% Question 2.3 % Dynamic Controller + LQR
h = 0.175;
t = 0;

Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
Gamma = [(exp(5*h)/16 + 5*exp(-3*h)/48-1/6); (1/3 -exp(-3*h)/3)];

Q = [0.01 0; 0 0.01]; R = 1;
Kb = dlqr(Phi, Gamma, Q, R);
Kb = [Kb 0.4]; % h = 0.175 (+ 0.1011 s)
% Kb = [Kb 0.35]; % h = 0.2 (+ 0.0878 s)
% Kb = [Kb 0.28]; % h = 0.25 (+ 0.0703 s)
% Kb = [Kb 0.2]; % h = 0.3 (+ 0.0502 s)
% Kb = [Kb 0.17]; % h = 0.35 (+ 0.0428 s)
tvec = 0:0.0001:h;
maxeig = zeros(3,length(tvec));
for i = 1:length(tvec)
    t = tvec(i);
    
    Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
    Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
    Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
   
    Phie = [[Phi, Gamma1]; zeros(1,3)];
    Gammae = [Gamma0; 1];
    
    maxeig(1,i) = max(abs(eig(Phie-Gammae*[K 0])));
    maxeig(2,i) = max(abs(eig(Phie-Gammae*Kb)));
end

plot(tvec, maxeig(1,:)); hold on; plot(tvec, maxeig(2,:)); plot(tvec, ones(1,length(tvec))); hold off
xlim([0,0.2]); grid on; grid minor;

tidx = find(maxeig(1,:) > 0.999 & maxeig(1,:) < 1.001);
maxK = tvec(round(mean(tidx)));

tidx = find(maxeig(2,:) > 0.999 & maxeig(2,:) < 1.001);
maxKb = tvec(round(mean(tidx)));
maxKb - maxK

%% Question 3.1
good = []; bad = [];
hvec = 0.01:0.01:0.2;
for i = 1:length(hvec)
    h = hvec(i);
    tvec = 0:0.001:h;
    for j = 1:length(tvec)
        t = tvec(j);
        if t+h <= 0.125
            Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
            Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
            Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
            Phie = [Phi, -Gamma0*K, -Gamma1*K;
                eye(size(Phi)),zeros(size(Phi)),zeros(size(Phi));
                zeros(size(Phi)),eye(size(Phi)),zeros(size(Phi))];
            L = eig(Phie); maxeig(i,j) = max(abs(L));
            %           Do this exactly like it's done in the lectures
            %             Phie = [[Phi, Gamma1, Gamma0]; 0 0 0 1; 0 0 0 0];
            %             Gammae = [0; 0; 0; 1];
            %             L = eig(Phie - Gammae*[K 0 0]); maxeig(i,j) = max(abs(L));
            if maxeig(i,j) < 1
                good = [good, [h;t+h]];
            else
                bad = [bad, [h;t+h]];
            end
        end
    end
end
scatter(good(1,:), good(2,:), 'bx');
scatter(bad(1,:), bad(2,:), 'rx'); hold off;

%% Question 3.3 % Dynamic Controller + LQR
h = 0.1; % h \in [0.07 , 0.1]
t = h;

Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
Gamma = [(exp(5*h)/16 + 5*exp(-3*h)/48-1/6); (1/3 -exp(-3*h)/3)];

Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];

Phie = [[Phi, Gamma1, Gamma0]; 0 0 0 1; 0 0 0 0];
Gammae = [0; 0; 0; 1];

Q = [0.01 0; 0 0.01]; R = 1;
Kb = dlqr(Phi, Gamma, Q, R);
% Kb = [Kb 0.15 0.15]; % h = 0.07 (+ 0.0164 s (max))
% Kb = [Kb 0.275 0.27]; % h = 0.1 (+ 0.0894 s (max))
Kb = [K 0.18 0.2]; % h = 0.1 (+ 0.0894 s (max))

tvec = 0:0.0001:h;
maxeig = zeros(2,length(tvec));
for i = 1:length(tvec)
    t = tvec(i);
    
    Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
    Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
    Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
    
    Phie = [[Phi, Gamma1, Gamma0]; 0 0 0 1; 0 0 0 0];
    Gammae = [0; 0; 0; 1];
    
    maxeig(1,i) = max(abs(eig(Phie-Gammae*[K 0 0])));
    maxeig(2,i) = max(abs(eig(Phie-Gammae*Kb)));
end

plot(tvec+h, maxeig(1,:)); hold on; plot(tvec+h, maxeig(2,:)); plot(tvec+h, ones(1,length(tvec))); hold off
xlim([h,2*h]); grid on; grid minor;

tidx = find(maxeig(1,:) > 0.999 & maxeig(1,:) < 1.001);
maxK = tvec(round(mean(tidx)));

tidx = find(maxeig(2,:) > 0.999 & maxeig(2,:) < 1.001);
if isnan(round(mean(tidx)))
    tf = length(tvec);
else
    tf = round(mean(tidx));
end
maxKb = tvec(tf);
maxKb - maxK

%% Question 4.1
K = [22.4 7];
h = 0.05;
tvec = [0.2*h; 0.5*h; h; 1.5*h];
clear CL
for i = 1:length(tvec)
    t = tvec(i);
    if t >= h
        t = t-h; % This structure is used for t = [h,2h), so we remove the offset.
    end
    Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
    Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
    Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
    if t < h
        CL{i} = [Phi-Gamma0*K, [0;0], Gamma1; 0, 0, 0, 1; -K, 0, 0];
    else
        CL{i} = [Phi, Gamma0, Gamma1; -K, 0, 0; 0, 0, 1, 0];
    end
end

n = 4;
cvx_begin sdp
variable P(n,n) symmetric
variable Q(n,n) symmetric
CL{1}'*P*CL{1} - P <= -Q;
CL{2}'*P*CL{2} - P <= -Q;
CL{3}'*P*CL{3} - P <= -Q;
CL{4}'*P*CL{4} - P <= -Q;
Q >= 10^(-3)*eye(n);
% P >= 10^(-2)*eye(n);
cvx_end
isposdef = [all(eig(P) > 0) all(eig(Q) > 0)]

%% Question 4.2
hvec = 0.001:0.001:0.35;
tvec = [0.2*h; 0.5*h; h; 1.5*h];
Ku = [0.15 0.15];
maxeig = zeros(2,length(tvec));
for i = 1:length(hvec)
    h = hvec(i);
    
    tt = 0.5*h
    if tt >= h
        t = tt-h;
    else
        t = tt;
    end
    Phi = [exp(5*h) 5/16*(exp(5*h)-exp(-3*h)); 0   exp(-3*h)];
    Gamma0 = [(5*exp(3*t - 3*h))/48 + exp(5*h - 5*t)/16 - 1/6; 1/3 - exp(3*t - 3*h)/3];
    Gamma1 = [(5*exp(-3*h))/48 + exp(5*h)/16 - (5*exp(-3*h)*exp(3*t))/48 - (exp(5*h)*exp(-5*t))/16; (exp(-3*h)*(exp(3*t) - 1))/3];
    
    if tt >= h
        CL = [Phi, Gamma1, Gamma0; 0, 0, 0, 1; -K, 0 0];   
        maxeig(1,i) = max(abs(eig(CL)));
        CL = [Phi, Gamma1, Gamma0; 0, 0, 0, 1; -K -Ku];
        maxeig(2,i) = max(abs(eig(CL)));    
    else
        CL = [Phi-Gamma0*K, [0;0], Gamma1; 0, 0, 0, 1; -K, 0, 0];
        maxeig(1,i) = max(abs(eig(CL)));
        CL = [Phi-Gamma0*K, [0;0], Gamma1-Gamma0*Ku(1); 0, 0, 0, 1; -K, 0, -Ku(1)];
        maxeig(2,i) = max(abs(eig(CL)));    
    end
    
end

tidx = find(maxeig(1,:) < 1.001);
maxK = hvec(max(tidx))

tidx = find(maxeig(2,:) < 1.001);
maxKb = hvec(max(tidx))
maxKb - maxK

plot(hvec, maxeig(1,:)); hold on; plot(hvec, maxeig(2,:)); plot(hvec, ones(1,length(hvec))); hold off
xlim([0,h]); grid on; grid minor;