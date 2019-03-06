function [tspan, X, U, Z, Y] = sim_kf_nlin(X0, Tp, Tend, AD, BD, C, d, Rf, Qf, K, Fd0, Fh0, Fc0, Tc, Th, Td, c, a, T0, h0)

P = eye(2);

tspan = 1:Tp:Tend;

x = X0;
z = [0;0];
X = zeros(length(tspan),2);
U = zeros(length(tspan),2);
Z = zeros(length(tspan),2);
Y = zeros(length(tspan),2);
for k=1:length(tspan)
    u = -K*z;
    [~,y] = ode15s(@(t,y)tank_derivatives_u_hold(y, u, d(k), Fd0, Fh0, Fc0, Tc, Th, Td, c, a),0:0.1:Tp, x + [h0;T0]);
    x = y(end,:)' - [h0;T0];
    ym = x + 1*randn(2,1);
    % predict
    z = AD*z + BD*u;
    P = AD*P*AD' + Qf;
    % update
    Kf = P*C'/(C*P*C' + Rf);
    z = z + Kf*(ym - z);
    P = P - Kf*C*P;
    X(k,:) = x;
    U(k,:) = u;
    Z(k,:) = z;
    Y(k,:) = ym;
end