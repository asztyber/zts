clear;
const_data
tspan = 0:.01:10;
X0 = [0; 0; 0.2; 0];
%u = 0;

% A = [0 1 0 0;
%     0 -b/M -m*g/M 0;
%     0 0 0 1;
%     0 b/(M*l) (m+M)*g/(M*l) 0];
% 
% B = [0; 1/M; 0; -1/(M*l)];

A = [0 1 0 0;
    0 -b*(I + m*l^2)/(I*(M+m) + m*M*l^2) -m^2*g*l^2/(I*(M+m) + m*M*l^2) 0;
    0 0 0 1;
    0 b*m*l/(I*(M+m) + m*M*l^2) (m+M)*g*m*l/(I*(M+m) + m*M*l^2) 0];

B = [0; (I + m*l^2)/(I*(M+m) + m*M*l^2); 0; -m*l/(I*(M+m) + m*M*l^2)];

Q = [1 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100];
R = 1;

K = lqr(A,B,Q,R);

%p = [-160 -8 -7 -9];%ok
%p = [-200 -18 -17 -19];%too agressive


%p = [-5 -1 -1.2 -1.5];%slow
%p = [-3 -1 -1.2 -1.5];%too slow

%K = place(A,B,p);

%% simulation
[t,X] = ode15s(@(t,X)pendulum_on_car_derivatives(X,-K*(X-[0; 0; 0; 0])),tspan,X0);

%plot(t,X)
%legend({'$\theta$', '$\dot{\theta}$', '$x$', '$\dot{x}$'}, 'Interpreter','latex')

xmin =  -1;
xmax = 5;
ymin =  -1;
ymax = 1;
figure(1);
axis([xmin xmax ymin ymax])
pbaspect([2.5 1 1])
hold;
hplot = plot(NaN,NaN,'-');
hrec = rectangle('Position',[X0(1)-w/2,-h/2,w,h]);
hcirc = rectangle('Position',[X0(1) + l*sin(X0(3)) - r, l*cos(X0(3)) - r, 2*r, 2*r],'Curvature',[1 1]);




for i = 1:length(X)
    plot_pendulum(X(i,:), w, h, l, r, hplot, hrec, hcirc)
    F(i) = getframe;
end

F(end+1) = getframe;
movie(F)
    