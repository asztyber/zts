%% const
clear;
c = 1.1;
a = 20;
Tc = 21;
Th = 77;
Td = 45;
Fc0 = 60;
Fh0 = 24;
Fd0 = 11;
tc = 100;
th = 180;
h0 = 22.56;
T0 = 37.93;
%% System parameters
Tp = 1;

A = [-2/(3*c*h0^3)*(Fh0 + Fc0 + Fd0) + a/(2*c*h0^(2.5)) 0;
    -3/(c*h0^4)*(Fh0*(Th - T0) + Fc0*(Tc - T0) + Fd0*(Td - T0)) -1/(c*h0^3)*(Fh0 + Fc0 + Fd0)];

B = [1/(3*c*h0^2) 1/(3*c*h0^2);
    1/(c*h0^3)*(Th - T0) 1/(c*h0^3)*(Tc - T0)];

Bd = [1/(3*c*h0^2); 1/(c*h0^3)*(Td - T0)];

C = eye(2);

D = zeros(2,3);
%% LQR controler
Q = eye(2);
R = eye(2);
[K,S,e] = lqr(A,B,Q,R);

%% transfer function
syms s;
G = K/(eye(2)*s - A)*B;
[nf, df] = numden(G(1,1));
tfn = sym2poly(nf);
tfd = sym2poly(df);
H = tf(tfn, tfd);
%nyquist(H)
% xlim([-.02 0.4])
% ylim([-.5 .5])

%% Kalman filter
Rf = B'*eye(2)*B;
Qf = 10^(-11)*eye(2);
L = lqr(A',C',Qf,Rf);
L = L';

%% LQG transfer function
syms s;
G2 = K/(eye(2)*s - A + B*K + L*C)*L*C/(eye(2)*s - A)*B;
[nf2, df2] = numden(G2(1,1));
tfn2 = sym2poly(nf2);
tfd2 = sym2poly(df2);
H2 = tf(tfn2, tfd2);
nyquist(H, H2)
xlim([-.2 0.2])
%ylim([-1.5 1.5])
legend('LQR','LQG', 'Location', 'northwest')