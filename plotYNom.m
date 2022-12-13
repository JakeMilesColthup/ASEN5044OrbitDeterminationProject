% Plot NL Measurements
figure()
sgtitle('Full Nominal Dynamics Simulation Measurment Data')

% Rho
subplot(4,1,1)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNom((ii*3) - 2,:), 'x')
end; clear ii;
ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Rhodot
subplot(4,1,2)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNom((ii*3) - 1,:))
end; clear ii;
ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% phi
subplot(4,1,3)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNom((ii*3),:), 's')
end; clear ii;
ylabel('$\phi ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% station IDs
subplot(4,1,4)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec, stationIDsMat(ii,:) , '^')
end; clear ii;
ylabel('Visible Station IDs')
xlabel('Time [s]')
