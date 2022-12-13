function xdot = nonLinearOrbitNoise(t,x,mu,wk)

r = sqrt(x(1)^2 + x(3)^2);
xdot(1) = x(2);
xdot(2) = (-mu * x(1) / r^3) + wk(2);
xdot(3) = x(4);
xdot(4) = (-mu * x(3) / r^3) + wk(4);
xdot = xdot';



end