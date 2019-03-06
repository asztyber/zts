function dX = tank_derivatives_u_hold(x, u, d, Fd0, Fh0, Fc0, Tc, Th, Td, c, a)

h = x(1);
T = x(2);

Fh = u(1) + Fh0;
Fc = u(2) + Fc0;
Fd = d + Fd0;


dX(1,1) = (Fh + Fc + Fd)/(3*c*h^2) -a/(3*c)*h^(-1.5);
dX(2,1) = (Fh*(Th - T) + Fc*(Tc - T) + Fd*(Td - T))/(c*h^3);