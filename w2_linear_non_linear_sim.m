clear;
Tend = 20;
tspan = 0:.1:Tend;
X0 = [0.1; 0]; % warunki poczatkowe
u = 0;
d = 0.1;

A = [0 1; -1 -d];
B = [0;1];
C = [1 0];
D = [0];

[t,X] = ode45(@(t,X)pendulum_derivatives(X,u,d),tspan,X0);

sys = ss(A, B, C, D);
resp = initial(sys,X0,tspan);

plot(t,X(:,1))
hold on;
plot(tspan,resp,'o')

title('$\theta(0) = 0.1$', 'Interpreter','latex')
legend({'$\theta$', '$\theta_d$ - model liniowy'}, 'Interpreter','latex')
