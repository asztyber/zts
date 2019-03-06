%% parameters
clear;
const_data
X0 = [0; 0; 0.1; 0];
Tp = 0.005;
Tend = 10;

A = [0 1 0 0;
    0 -b*(I + m*l^2)/(I*(M+m) + m*M*l^2) -m^2*g*l^2/(I*(M+m) + m*M*l^2) 0;
    0 0 0 1;
    0 b*m*l/(I*(M+m) + m*M*l^2) (m+M)*g*m*l/(I*(M+m) + m*M*l^2) 0];

B = [0; (I + m*l^2)/(I*(M+m) + m*M*l^2); 0; -m*l/(I*(M+m) + m*M*l^2)];

C = [1 0 0 0;
    0 0 1 0];

%C = [1 0 0 0];

Q = [1 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100];
R = 0.001;

%% LQR controler
[K,S,e] = lqr(A,B,Q,R);

%% transfer function
syms s;
G = K/(eye(4)*s - A)*B;
[nf, df] = numden(G);
tfn = sym2poly(nf);
tfd = sym2poly(df);
H = tf(tfn, tfd);
nyquist(H)

%% Kalman filter
Qf = 10^(20)*eye(4);
Rf = eye(2);
L = lqr(A',C',Qf,Rf);
L = L';

%% LQG transfer function
syms s;
G = K/(eye(4)*s - A + B*K + L*C)*L*C/(eye(4)*s - A)*B;
[nf, df] = numden(G);
tfn = sym2poly(nf);
tfd = sym2poly(df);
H = tf(tfn, tfd);
nyquist(H)
%xlim([-1.2 0.05])
%ylim([-1.5 1.5])