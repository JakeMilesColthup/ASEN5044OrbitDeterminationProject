% Plot Noisy Simulated Ground Truth States
figure()
sgtitle('Typical Simulation LKF: Noisy Ground Truth States')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xNL(1,:),'LineWidth',1.5)
plot(tvec,xTMT(1,:, 1),'LineWidth',1.5)
plot(tvec,xLKF(1,:, 1),'LineWidth',1.5)
plot(tvec,xEKF(1,:, 1),'LineWidth',1.5) 
ylabel('X [km]')
xlabel('Time [s]')
legend('NL','TMT','LKF','EKF')

% Xdot
subplot(4,1,2) 
grid on; box on; hold on;
plot(tvec,xNL(2,:),'LineWidth',1.5)
plot(tvec,xTMT(2,:, 1),'LineWidth',1.5)
plot(tvec,xLKF(2,:, 1),'LineWidth',1.5)
plot(tvec,xEKF(2,:, 1),'LineWidth',1.5) 
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xNL(3,:),'LineWidth',1.5)
plot(tvec,xTMT(3,:, 1),'LineWidth',1.5)
plot(tvec,xLKF(3,:, 1),'LineWidth',1.5)
plot(tvec,xEKF(3,:, 1),'LineWidth',1.5) 
ylabel('Y [km]')
xlabel('Time [s]')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xNL(4,:),'LineWidth',1.5)
plot(tvec,xTMT(4,:, 1),'LineWidth',1.5)
plot(tvec,xLKF(4,:, 1),'LineWidth',1.5)
plot(tvec,xEKF(4,:, 1),'LineWidth',1.5) 
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
