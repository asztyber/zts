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

%% Tustin discretization
AD = (eye(2) + 1/2*A*Tp)/((eye(2) - 1/2*A*Tp));
BD = A\(AD - eye(2))*B;
BdD = A\(AD - eye(2))*Bd;

%% LQR discrete time, infinite horizon
Q = eye(2);
R = eye(2);
R1 = 0.05*eye(2);

[K,S,e] = lqrd(A,B,Q,R,Tp); % continuous time matrices
[K1,S1,e1] = lqrd(A,B,Q,R1,Tp);
eigs(AD - BD*K1)
%% LQR discrete time, finite horizon
Tend = 2000;
P = lqr_p_list(AD,BD,Q,R1,Q,Tend);
Pflat = reshape(P,4,Tend);
subplot(2,2,1)
plot(Pflat(1,:))
hold on
plot(S1(1,1)*ones(Tend,1))
legend({'$P_\tau^{11}$', '$P_\infty^{11}$'}, 'Interpreter','latex','Location','southwest' )
subplot(2,2,2)
plot(Pflat(2,:))
hold on
plot(S1(1,2)*ones(Tend,1))
legend({'$P_\tau^{21}$', '$P_\infty^{21}$'}, 'Interpreter','latex','Location','northwest' )
subplot(2,2,3)
plot(Pflat(3,:))
hold on
plot(S1(2,1)*ones(Tend,1))
legend({'$P_\tau^{12}$', '$P_\infty^{12}$'}, 'Interpreter','latex','Location','northwest' )
subplot(2,2,4)
plot(Pflat(4,:))
hold on
plot(S1(2,2)*ones(Tend,1))
legend({'$P_\tau^{22}$', '$P_\infty^{22}$'}, 'Interpreter','latex','Location','southwest' )
%% simulation
%d=cumsum(0.2*randn(Tend,1));
d=10*ones(Tend,1);
X0 = [-10;5];

%[tspan, X, U] = sim_state_feedback(X0, Tp, Tend, AD, BD, BdD, d, K);
%[tspan, X1, U1] = sim_state_feedback(X0, Tp, Tend, AD, BD, BdD, d, K1);
%[tspan, X, U] = sim_state_feedback_nlin(X0, Tp, Tend, d, zeros(2,2), Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0);
[tspan, X1, U1] = sim_state_feedback_nlin(X0, Tp, Tend, d, K1, Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0);
[tspan, Xk, Uk] = sim_p_nlin(X0, Tp, Tend, AD, BD, d, P, R1, Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0);
%% plots x
plot(X1(:,1), 'b')
hold on
plot(X1(:,2), 'r')
hold on
plot(Xk(:,1), 'b--')
hold on
plot(Xk(:,2), 'r--')
legend({'$h$, $P_\infty$', '$T$, $P_\infty$','$h$, $P_k$', '$T$, $P_k$'}, 'Interpreter','latex','Location','southeast' )
%% plots u
plot(U1(:,1), 'r')
hold on
plot(U1(:,2), 'b')
hold on
plot(Uk(:,1), 'r--')
hold on
plot(Uk(:,2), 'b--')
legend({'$F_H$, $P_\infty$', '$F_C$, $P_\infty$','$F_H$, $P_k$', '$F_C$, $P_k$'}, 'Interpreter','latex','Location','northeast' )