function plot_pendulum(X, w, h, L, r, hplot, hrec, hcirc)

theta = X(3);
x = X(1);

set(hplot,'XData',[x x + L*sin(theta)],'YData',[0 L*cos(theta)]);
set(hrec, 'Position', [x-w/2,-h/2,w,h])
set(hcirc, 'Position',[x + L*sin(theta) - r, L*cos(theta) - r, 2*r, 2*r])

