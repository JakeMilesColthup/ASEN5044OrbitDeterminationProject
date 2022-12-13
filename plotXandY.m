% Plot Noisy Simulated Ground Truth States
figure()
sgtitle('Typical Simulation LKF: Noisy Ground Truth X and Y')

% X
grid on; box on; hold on;
plot(xNL(1,1:nTimesteps), xNL(3,1:nTimesteps),'LineWidth',1.5)
plot(xTMT(1,1:nTimesteps, 1), xTMT(3,1:nTimesteps, 1), 'LineWidth',1.5)
ylabel('Y [km]')
xlabel('X [km]')
legend('NL','TMT')
