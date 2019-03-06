function K = kalman_filter(Ad,C,Qf,Rf)
% initial values
P = eye(4);
K = P*C'/(C*P*C' + Rf);

for i=1:100
    P = Ad*P*Ad' + Qf;
    K = P*C'/(C*P*C' + Rf);
    P = P - K*C*P;
end
