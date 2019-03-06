function [tspan, X, U] = sim_state_feedback_lin(X0, Tp, Tend, AD, BD, BdD, d, K)
tspan = 1:Tp:Tend;

x = X0;
X = zeros(length(tspan),2);
U = zeros(length(tspan),2);
for k=1:length(tspan)
    u = -K*x;
    x = AD*x + BD*u + BdD*d(k);
    [tn,yn] = ode15s(@(t,x)tank_derivatives(x, t, u, Fd0, Fh0, Fc0, Tc, Th, Td, c, a),t,X0);
    X(k,:) = x;
    U(k,:) = u;
end