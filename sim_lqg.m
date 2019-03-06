function [tspan, X, U, Z, Y] = sim_lqg(X0, Tp, Tend, K, Ad, Bd, C, L)
tspan = 0:Tp:Tend;

x = X0;
%z = X0;
z = [0;0;0;0];
X = zeros(length(tspan),4);
U = zeros(length(tspan),1);
Z = zeros(length(tspan),4);
Y = zeros(length(tspan),4);
for k=1:length(tspan)
    u = -K*z;
    %u = -K*x;
    [~,y] = ode15s(@(t,y)pendulum_on_car_derivatives(y,u),0:Tp/10:Tp,x);
    x = y(end,:)';
    y = x + 0.01*randn(4,1);
    %y = x;
    z = Ad*z + Bd*u + L*(C*y - C*z);
    X(k,:) = x;
    U(k,:) = u;
    Z(k,:) = z;
    Y(k,:) = y;
end