%% parameters
clear;
const_data
X0 = [0; 0; 0.2; 0];
Tp = 0.005;
Tend = 10;

A = [0 1 0 0;
    0 -b*(I + m*l^2)/(I*(M+m) + m*M*l^2) -m^2*g*l^2/(I*(M+m) + m*M*l^2) 0;
    0 0 0 1;
    0 b*m*l/(I*(M+m) + m*M*l^2) (m+M)*g*m*l/(I*(M+m) + m*M*l^2) 0];

B = [0; (I + m*l^2)/(I*(M+m) + m*M*l^2); 0; -m*l/(I*(M+m) + m*M*l^2)];

%C = [1 0 0 0;
%    0 0 1 0];

C = [1 0 0 0];

Q = [1 0 0 0;
    0 1 0 0;
    0 0 10 0;
    0 0 0 100];
R = 0.001;

%% exact
Ad = expm(A*Tp); % exp - elementwise
Bd = A\(expm(A*Tp) - eye(4))*B;
Bd = [0;Bd(2:4)];
%% LQR controler
[K,S,e] = lqrd(A,B,Q,R,Tp);

%% observability
obsv(A,C)
rank(obsv(A,C))

%% observer
%L_T=place(Ad',C',[0.3,0.1,0.15,0.05]);
%L_T=place(Ad',C',[0.8,0.9,0.85,0.87]);
L_T=place(Ad',C',[0.9,0.95,0.98,0.97]);
L=L_T';

%% Kalman Filter
%Qf = 10^(-4)*eye(4);
Qf = [10^(-6) 0 0 0;
    0 10^(-6) 0 0;
    0 0 10^(-4) 10^(-4);
    0 0 10^(-4) 10^(-2)];
Rf = 0.001*eye(2);
[P,Lf,Gf] = dare(Ad',C',Qf,Rf);
Kf = P*C'/(C*P*C' + Rf);
%Kf = kalman_filter(Ad,C,Qf,Rf);
eigs(Ad - Kf*C)
%% simulation
[tspan, X, U, Z, Y] = sim_lqg(X0, Tp, Tend, K, Ad, Bd, C, L);

%% plots
axlim=2000;
subplot(2,2,1)
plot(X(:,1));
hold on 
plot(Z(:,1),'g')
plot(Y(:,1),'.');
xlim([0 axlim])
legend({'$x$','$\hat{x}$', '$x_M$',},'Interpreter','latex', 'Location','southeast')
subplot(2,2,2)
plot(X(:,2))
hold on 
plot(Z(:,2),'g')
xlim([0 axlim])
legend({'$\dot{x}$', '$\hat{\dot{x}}$'},'Interpreter','latex')
subplot(2,2,3)
plot(X(:,3))
hold on 
plot(Z(:,3),'g')
plot(Y(:,3),'.')
xlim([0 axlim])
legend({'$\theta$', '$\hat{\theta}$', '$\theta_M$'},'Interpreter','latex', 'Location','northeast')
subplot(2,2,4)
plot(X(:,4))
hold on 
plot(Z(:,4),'g')
xlim([0 axlim])
legend({'$\dot{\theta}$', '$\hat{\dot{\theta}}$'},'Interpreter','latex', 'Location','northeast')

%% picture

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