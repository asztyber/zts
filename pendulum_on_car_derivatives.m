function dX = pendulum_on_car_derivatives(X,u)

const_data

theta = X(3);
dtheta = X(4);
dx = X(2);

D = (m+M)*(m*l^2 + I) - m^2*l^2*cos(theta)^2;
ddx = u*(m*l^2 + I)/D + ((m*l*sin(theta)*dtheta^2 - b*dx)*(m*l^2 + I) - g*m^2*l^2*cos(theta)*sin(theta))/D;
ddtheta = m*l*cos(theta)/D*(-u - m*l*sin(theta)*dtheta*2 + b*dx) + m*g*l*sin(theta)/(m*l^2 + I) + g*m^2*l^2*cos(theta)*sin(theta)/D/(m*l^2 + I);

dX(1,1) = dx;
dX(2,1) = ddx;
dX(3,1) = dtheta;
dX(4,1) = ddtheta;