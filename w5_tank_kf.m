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
Tp = 5;

A = [-2/(3*c*h0^3)*(Fh0 + Fc0 + Fd0) + a/(2*c*h0^(2.5)) 0;
    -3/(c*h0^4)*(Fh0*(Th - T0) + Fc0*(Tc - T0) + Fd0*(Td - T0)) -1/(c*h0^3)*(Fh0 + Fc0 + Fd0)];

B = [1/(3*c*h0^2) 1/(3*c*h0^2);
    1/(c*h0^3)*(Th - T0) 1/(c*h0^3)*(Tc - T0)];

Bd = [1/(3*c*h0^2); 1/(c*h0^3)*(Td - T0)];

C = eye(2);

D = zeros(2,3);

%% Tustin discretization
AD = (eye(2) + 1/2*A*Tp)/((eye(2) - 1/2*A*Tp));
BD = A\(AD - eye(2))*B;
BdD = A\(AD - eye(2))*Bd;

%% LQR discrete time, infinite horizon
Q = eye(2);
R = eye(2);

[K,S,e] = lqrd(A,B,Q,R,Tp); % continuous time matrices
eigs(AD - BD*K);

%% Kalman Filter
Rf = eye(2);
Qf = Bd'*eye(2)*Bd;

%% simulation
Tend = 2000;
%d=cumsum(1*randn(Tend,1));
%d=10*ones(Tend,1);
d = randn(Tend,1);
X0 = [-10;5];

[tspan, X, U, Z, Y] = sim_kf_nlin(X0, Tp, Tend, AD, BD, C, d, Rf, Qf, K, Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0);

%% plots
plot(X(:,1), 'b')
hold on
plot(X(:,2), 'r')
plot(Z(:,1), 'b--')
plot(Z(:,2), 'r--')
plot(Y(:,1), 'b.')
plot(Y(:,2), 'r.')
legend({'$h$', '$T$','$\hat{h}$', '$\hat{T}$', '$h_M$', '$T_M$'}, 'Interpreter','latex','Location','southeast' )