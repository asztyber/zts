function dX = pendulum_derivatives(X,u,d)

theta = X(1);
dtheta = X(2);

dX(1,1) = dtheta;
dX(2,1) = -sin(theta) - d*dtheta;
