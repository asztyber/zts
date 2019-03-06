function dX = tank_derivatives(x, t, u, Fd0, Fh0, Fc0, Tc, Th, Td, c, a)

h = x(1);
T = x(2);

if t < 1
    Fh = Fh0;
    Fc = Fc0;
    Fd = Fd0;
else
    Fh = interp1(u(:,1),t) + Fh0;
    Fc = interp1(u(:,2),t) + Fc0;
    Fd = interp1(u(:,3),t) + Fd0;
end

dX(1,1) = (Fh + Fc + Fd)/(3*c*h^2) -a/(3*c)*h^(-1.5);
dX(2,1) = (Fh*(Th - T) + Fc*(Tc - T) + Fd*(Td - T))/(c*h^3);

% V = x(1);
% h = nthroot(V/c,3);
% T = x(2);
% Fh = u(1);
% Fc = u(2);
% 
% dX(1,1) = Fh + Fc + Fd -a*sqrt(h);
% dX(2,1) = (Fh*(Th - T) + Fc*(Tc - T) + Fd*(Td - T))/(c*h^3);