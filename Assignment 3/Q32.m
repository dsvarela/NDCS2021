close all; clearvars;
H = .00100*rand(4); H = H*H';

H11 = H(1:2,1:2);
H12 = H(1:2,2:4);
H12 = H(1:2,3:4);
H21 = H(3:4,1:2);
H22 = H(3:4,3:4);

syms lam
figure('Position', [10, 50, 950,900])
fplot(det(lam * eye(2) - H11\ H12 /H22* H21/lam)); 
hold on; viscircles([0,0],1);
yline(0, 'r', 'linewidth',2)
xlim([-2,2]); ylim([-2,2]);
grid on; grid minor;
hold off;

figure('Position', [900, 50, 950,900])
fplot(det(lam * H11 -  H12 /H22* H21/lam)); 
hold on; viscircles([0,0],1);
yline(0, 'r', 'linewidth',2)
xlim([-2,2]); ylim([-1e-6,1e-6]);
grid on; grid minor;
hold off;


A = [zeros(2,2),H11\H12; H22\H21, zeros(2,2)];
[eig(H), eig(A)]; [ eig(H11), eig(H12 /H22* H21)];
eig(eye(2) - H11\ H12 /H22* H21)
    eig(H11 - H12 /H22* H21)
