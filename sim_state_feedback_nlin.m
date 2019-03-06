function [tspan, X, U] = sim_state_feedback_nlin(X0, Tp, Tend, d, K, Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0)
tspan = 1:Tp:Tend;

x = X0;
X = zeros(length(tspan),2);
U = zeros(length(tspan),2);
for k=1:length(tspan)
    u = -K*x;
    [~,y] = ode15s(@(t,y)tank_derivatives_u_hold(y, u, d(k), Fd0, Fh0, Fc0, Tc, Th, Td, c, a),0:0.1:Tp, x + [h0;T0]);
    x = y(end,:)' - [h0;T0];
    X(k,:) = x;
    U(k,:) = u;
end