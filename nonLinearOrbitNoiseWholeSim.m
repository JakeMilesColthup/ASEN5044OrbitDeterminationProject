function xdot = nonLinearOrbitNoiseWholeSim(t,x,mu,wk, tvec)

r = sqrt(x(1)^2 + x(3)^2);
xdot(1) = x(2);
xdot(2) = (-mu * x(1) / r^3);
xdot(3) = x(4);
xdot(4) = (-mu * x(3) / r^3);
xdot = xdot';

diffVec = abs(t-tvec);
log = diffVec < 1e-1;

if any(log)
    xdot = xdot + wk(:, log);
end % if

end