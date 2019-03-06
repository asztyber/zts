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

V0 = h0^3*c;
X0 = [h0 T0];

A = [-2/(3*c*h0^3)*(Fh0 + Fc0 + Fd0) + a/(2*c*h0^(2.5)) 0;
    -3/(c*h0^4)*(Fh0*(Th - T0) + Fc0*(Tc - T0) + Fd0*(Td - T0)) -1/(c*h0^3)*(Fh0 + Fc0 + Fd0)];

B = [1/(3*c*h0^2) 1/(3*c*h0^2);
    1/(c*h0^3)*(Th - T0) 1/(c*h0^3)*(Tc - T0)];

Bd = [1/(3*c*h0^2); 1/(c*h0^3)*(Td - T0)];

C = eye(2);

D = zeros(2,3);

sys = ss(A, [B Bd], C, D, 'InputDelay', [0 0 0]);

n = 8000;
t= 0:n;
Fhin = [-10*ones(round(n/3),1);zeros(round(n/3),1);ones(round(n/3),1)];
Fdin = cumsum(0.05*randn(n+1,1));
u = [Fhin zeros(n+1,1) Fdin];
ylin =lsim(sys,u,t);

[tn,yn] = ode15s(@(t,x)tank_derivatives(x, t, u, Fd0, Fh0, Fc0, Tc, Th, Td, c, a),t,X0);

subplot(2,1,1)
plot(Fhin)
hold on
plot(Fdin)
xlim([0,8000])
legend({'$F_{H} - F_{H0}$', '$F_D - F_{D0}$',}, 'Interpreter','latex','Location','southeast' )
subplot(2,1,2)
plot(tn,yn(:,2))
hold on;
plot(t,ylin(:,2) + T0)
legend({'$T$', '$T$ - linearized',}, 'Interpreter','latex','Location','southeast' )