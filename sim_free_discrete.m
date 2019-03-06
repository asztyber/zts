function [tspand, Xd] = sim_free_discrete(X0, Tp, Tend, Ad)

xd = X0;
tspand = 0:Tp:Tend;
Xd = zeros(length(tspand),2);
for k=1:length(tspand)
    xd = Ad*xd; % + Bd*u
    Xd(k,:) = xd';
end