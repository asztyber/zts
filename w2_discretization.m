%% discretization
clear;
Tend = 10;
tspan = 0:.01:Tend;
X0 = [0.1; 0];
u = 0;
d = 0.1;
Tp = 0.1;

A = [0 1; -1 -d];
B = [0;1];
C = [1 0];
D = [0];

% exact
Ad = expm(A*Tp); % exp - elementwise
Bd = A\(expm(A*Tp) - eye(2))*B;

% Euler forward
 Ad_ef = eye(2) + A*Tp;
 Bd_ef = B*Tp;

% Euler backward
 Ad_eb = inv(eye(2) - A*Tp);
 Bd_eb = A\(Ad_eb - eye(2))*B;

% Tustin
Ad_t = (eye(2) + 1/2*A*Tp)/((eye(2) - 1/2*A*Tp));
Bd_t = A\(Ad_t - eye(2))*B;

% Tustin - Matlab library
sys = ss(A, B, C, D);
sysd_tustin = c2d(sys,Tp,'tustin');

%% simulation
[t,X] = ode45(@(t,X)pendulum_derivatives(X,u,d),tspan,X0);
[tspand, Xd_ef] = sim_free_discrete(X0, Tp, Tend, Ad_ef);
[tspand, Xd_eb] = sim_free_discrete(X0, Tp, Tend, Ad_eb);
[tspand, Xd_t] = sim_free_discrete(X0, Tp, Tend, Ad_t);
[tspand, Xd] = sim_free_discrete(X0, Tp, Tend, Ad);

yd_approx = [X0(1); X0(1); zeros(length(tspand),1)];
for k=3:(length(tspand)+2)
    yd_approx(k) = Tp^2*u + (2 - d*Tp)*yd_approx(k - 1) + (-1 + d*Tp - Tp^2)*yd_approx(k - 2);
end

plot(t,X(:,1))
hold on;
plot(tspand,Xd(:,1),'.')
plot(tspand,Xd_t(:,1),'--')
plot(tspand,Xd_eb(:,1),'--')
plot(tspand,Xd_ef(:,1),'o')
plot(tspand,yd_approx(2:end - 1))
legend({'$\theta$', 'exact', 'Tustin', 'Euler backward', 'Euler forward', 'Euler forward - z'}, 'Interpreter','latex')
title('$T_P = 0.01$', 'Interpreter','latex')