% ASEN 5044
% Billy Weigl, Jake Miles, Adam Buencamino
% Project

% Housekeeping
clear; clc; close all;

%% Part 1

%% Read in Data
% Load data
load('orbitdeterm_finalproj_KFdata.mat')

% Store number of timesteps
timesteps = length(tvec);

%% Problem 2/3 Parameter Initialization

% Note: since the nominal trajectory changes wrt time, the corresponding lineraized DT model
% is not LTI, it is time-varying.

% Therefore, the calculated DT model is updated and changes every single time, 
% but the derived generalized version of the matrices can be found in the report
% with the values of x not evaluated (since those change with time).

% Define parameters
mu = 398600; % km^3/s^2
deltaT = 10; % s
RE = 6378; % km
omegaE = (2*pi)/86400; % rad/s
nStations = 12;

% Define nominal operating point
r0 = 6678; % km
x0 = [6678; 0; 0; r0*sqrt(mu/r0^3)]; % km & km/s

% Define perturbation for part 1
% perturb = [0; 0.075; 0; -0.021];

% Define perturbation for the rest 
perturb = [.001 ; 0.0005; -.0004; -0.0002];

x0_pert = x0 + perturb;

% v0 = r0*sqrt(mu/r0^3); % km/s
% omega0 = v0/r0; % rad/s
% ustar = [0; 0]; % Nominal control inputs should be 0 


%% Problem 2/3 Computation

%% First, simulate full nonlinear dynamics using ode45 to get nominal trajectory

% Define time vector 
tspan = tvec(1:end);

% Define ode tolerance
tolPts = odeset('RelTol',5e-14,'AbsTol',1e-25);
% Prof used e-12 e-12

% Call ode45 to find nominal trajectory
[tOde,xNom] = ode45(@(t_total,x_total) nonLinearOrbit(t_total,x_total,mu), tspan, x0, tolPts);

% Rearrange output
xNom = interp1(tOde, xNom, tvec);
xNom = xNom';

% Call ode45 with perturbations to simulate the full nonlinear dynamics
[tOde,xNL] = ode45(@(t_total,x_total) nonLinearOrbit(t_total,x_total,mu), tspan, x0_pert, tolPts);

% Rearrange output
xNL = interp1(tOde, xNL, tvec);
xNL = xNL';


% Plot Full Nonlinear Simulation Results
figure()
sgtitle('Full Nonlinear Dynamics Simulation State History')

% X
subplot(4,1,1)
grid on; box on;
plot(tvec,xNL(1,:),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')

% Xdot
subplot(4,1,2)
grid on; box on;
plot(tvec,xNL(2,:),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on;
plot(tvec,xNL(3,:),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')

% Xdot
subplot(4,1,4)
grid on; box on;
plot(tvec,xNL(4,:),'LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')


%% Simulate linearized dynamics

% Initialize remaining states and predicted output vectors
xLinearized = zeros(length(x0_pert), length(tvec));
xLinearized(:,1) = x0_pert;
yLinearized = NaN * ones(3*nStations, length(tvec)); % NaN matrix, will only be overwritten if the station can see the satellite
yNL = NaN * ones(3*nStations, length(tvec));
yNom = NaN * ones(3*nStations, length(tvec));
stationIDsMat = NaN * ones(nStations, length(tvec));
Ftilde = zeros(length(x0_pert), length(x0_pert), length(tvec));
Omegatilde = zeros(length(x0_pert), 2, length(tvec));
Htilde = NaN * ones(3*nStations, length(x0_pert), length(tvec));

[Ftilde(:,:,1), ~, Omegatilde(:,:,1)] = computeLinearizedDyn(xNom(:,1), deltaT, mu);
[stationVec, stationIDs] = checkVisibleStations(nStations, RE, omegaE, xNL(:,1), tvec(1));
Htilde(:,:,1) = computeLinearizedH(xNom(:,1), stationVec);

% Simulate linearized DT dynamics and measurement models
for ii = 1:length(tvec)-1

    % Linearize system and find F, G, Omega
    [Ftilde(:,:,ii+1), Gtildek, Omegatilde(:,:,ii+1)] = computeLinearizedDyn(xNom(:,ii), deltaT, mu); 

    % Simulate DT Dynamics
    xLinearized(:,ii+1) = xNom(:,ii+1) + Ftilde(:,:,ii+1)*( xLinearized(:,ii) - xNom(:,ii) ) + Gtildek*[0;0] + Omegatilde(:,:,ii+1)*[0;0]; % Note: no control input perturbations, no process noise

    % Calculate information of the visible stations
    [stationVec, stationIDs] = checkVisibleStations(nStations, RE, omegaE, xNL(:,ii+1), tvec(ii+1));
    [stationVecNom, stationIDsNom] = checkVisibleStations(nStations, RE, omegaE, xNom(:,ii+1), tvec(ii+1));

    % Calculate H
    Htilde(:,:,ii+1) = computeLinearizedH(xNom(:,ii+1), stationVec);

    % Simulate measurement model
    yNL(:,ii+1) = computeYNL(nStations, stationIDs, RE, omegaE, xNL(:,ii+1), tvec(ii+1));
    yNom(:,ii+1) = computeYNL(nStations, 1:nStations, RE, omegaE, xNom(:,ii+1), tvec(ii+1));
    yLinearized(:,ii+1) = yNom(:,ii+1) + Htilde(:,:,ii+1) * ( xLinearized(:,ii+1) - xNom(:,ii+1) ); % No measurement noise

    % Store which station IDS are visible
    for jj = 1:length(stationIDs)
        stationIDsMat(stationIDs(jj), ii+1) = stationIDs(jj);
    end % for

end; clear ii; % for

% Determine perturbation from nominal trajectory at all times
xPerturbs = xLinearized - xNom; 


%% Plot Results

% Plot NL Measurements
figure()
sgtitle('Full Nonlinear Dynamics Simulation Measurment Data')

% Rho
subplot(4,1,1)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNL((ii*3) - 2,:), 'x')
end; clear ii;
ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Rhodot
subplot(4,1,2)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNL((ii*3) - 1,:))
end; clear ii;
ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% phi
subplot(4,1,3)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNL((ii*3),:), 's')
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


% Plot Linearized Dynamics Simulation Results
figure()
sgtitle('Linearized Dynamics Simulation State History')

% X
subplot(4,1,1)
grid on; box on;
plot(tvec,xLinearized(1,:),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')

% Xdot
subplot(4,1,2)
grid on; box on;
plot(tvec,xLinearized(2,:),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on;
plot(tvec,xLinearized(3,:),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')

% Ydot
subplot(4,1,4)
grid on; box on;
plot(tvec,xLinearized(4,:),'LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')


% Plot Both Sim Results on Same Plot
figure()
sgtitle('State History for Both Simulations')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xNL(1,:),'LineWidth',1.5)
plot(tvec,xLinearized(1,:),'--','LineWidth',2.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xNL(2,:),'LineWidth',1.5)
plot(tvec,xLinearized(2,:),'--','LineWidth',2.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xNL(3,:),'LineWidth',1.5)
plot(tvec,xLinearized(3,:),'--','LineWidth',2.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xNL(4,:),'LineWidth',1.5)
plot(tvec,xLinearized(4,:),'--','LineWidth',2.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')


% Plot Linearized Dynamics Simulation Perturbations
figure()
sgtitle('Linearized Dynamics Simulation Perturbations')

% X
subplot(4,1,1)
grid on; box on;
plot(tvec,xPerturbs(1,:),'LineWidth',1.5)
ylabel('$\delta{X}$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Xdot
subplot(4,1,2)
grid on; box on;
plot(tvec,xPerturbs(2,:),'LineWidth',1.5)
ylabel('$\delta{\dot{X}}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on;
plot(tvec,xPerturbs(3,:),'LineWidth',1.5)
ylabel('$\delta{Y}$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Ydot
subplot(4,1,4)
grid on; box on;
plot(tvec,xPerturbs(4,:),'LineWidth',1.5)
ylabel('$\delta{\dot{Y}}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')


% Plot Linearized Measurements
figure()
sgtitle('Linearized Dynamics Simulation Measurment Data')

% Rho
subplot(4,1,1)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yLinearized((ii*3) - 2,:), 'x')
end; clear ii;
ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Rhodot
subplot(4,1,2)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yLinearized((ii*3) - 1,:))
end; clear ii;
ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% phi
subplot(4,1,3)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yLinearized((ii*3),:), 's')
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


%% Part 2

%% Problem 4

%% Tune Initialization of LKF

% After extensive tuning, the following initial covariance matrix was selected
P0 = diag([1e-10 1e-10 1e-10 1e-10]);

Pplus0 = P0;
P0_noise = diag([0 Qtrue(1,1) 0 Qtrue(2,2)]);

% Set seed just for help debugging/analyzing results
rng(100);

% Attempt of initalizing delxplus0
delxplus0 = perturb;


%% Generate TMT Simulation Data
% Define Monte Carlo Parameters 
nRuns = 25; % Number of Monte Carlo runs
nTimesteps = 100; % 100 timesteps corresponds to 1000 s

% Initialize output truth model data
xTMT = zeros(4, nTimesteps, nRuns);

% Run Monte Carlo Sims and extract TMT data
for ii = 1:nRuns
    xTMT(:,:,ii) = simulateTruthModel(nTimesteps, Qtrue, x0_pert, P0_noise, deltaT);
end; clear ii; % for

% Adjust timevector 
tvec = tvec(1:nTimesteps);


%% Plot Noisy Simulated Ground Truth States for One Instance

% Plot Noisy Simulated Ground Truth States
figure()
sgtitle('Typical Simulation LKF: Noisy Ground Truth States')

% X
subplot(4,1,1)
grid on; box on;
plot(tvec,xTMT(1,:, 1),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')

% Xdot
subplot(4,1,2)
grid on; box on;
plot(tvec,xTMT(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on;
plot(tvec,xTMT(3,:, 1),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')

% Ydot
subplot(4,1,4)
grid on; box on;
plot(tvec,xTMT(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')


%% Simulate Noisy Measurements

% Initialize output  data
yNoisy = NaN * ones(3*nStations, nTimesteps, nRuns);

% Simulate noisy measurements with created function 
for ii = 1:nRuns
    yNoisy(:,:,ii) = simulateNoisyMeas(yNL(:,1:nTimesteps), Rtrue, nStations);
end; clear ii; % for


%% Plot Noisy Simulated Data

% Plot NL Measurements
figure()
sgtitle('Typical Simulation LKF: Noisy Measurements')

% Rho
subplot(4,1,1)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisy((ii*3) - 2,:, 1), 'x')
end; clear ii;
ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Rhodot
subplot(4,1,2)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisy((ii*3) - 1,:, 1))
end; clear ii;
ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% phi
subplot(4,1,3)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisy((ii*3),:, 1), 's')
end; clear ii;
ylabel('$\phi ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% station IDs
subplot(4,1,4)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec, stationIDsMat(ii,1:nTimesteps) , '^')
end; clear ii;
ylabel('Visible Station IDs')
xlabel('Time [s]')


%% Run the LKF

% Adjust sizes for consistency
xNom = xNom(:,1:nTimesteps);

% Initialize outputs
xPerturbsLKF = zeros(4, nTimesteps, nRuns);
xLKF = zeros(4, nTimesteps, nRuns);
xErrorLKF = zeros(4, nTimesteps, nRuns);
exLKF = zeros(4, nTimesteps, nRuns);
PLKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nRuns); % 4x4x1401xnRuns - Combined so not dependent on stations
SLKF = cell(nTimesteps, nRuns);
yErrorLKF = cell(nTimesteps, nRuns);
eyLKF = cell(nTimesteps, nRuns); % Perturbed innovation
delyLKF = cell(nTimesteps, nRuns); % Meas. innovation

% Run LKF
for ii = 1:nRuns
    [xPerturbsLKF(:,:,ii), xLKF(:,:,ii), xErrorLKF(:,:,ii), PLKF(:,:,:,ii), SLKF(:,ii), yErrorLKF(:,ii), exLKF(:,:,ii), eyLKF(:,ii), delyLKF(:,ii)] = runLKFCombined(xTMT(:,:,ii), xNom, yNom, yNoisy(:,:,ii), Qtrue, Rtrue, Ftilde, Omegatilde, Htilde, delxplus0, Pplus0);
end; clear ii; % for


%% Calculate 2sigma Bounds

% Initialize output
twoSigmasLKF = zeros(4, nTimesteps, nRuns);

% Loop through and calculate 2 sigma bounds
for jj = 1:nRuns
    for ii = 1:nTimesteps
        twoSigmasLKF(:, ii, jj) = [2*sqrt(PLKF(1,1,ii,jj)); 2*sqrt(PLKF(2,2,ii,jj)); 2*sqrt(PLKF(3,3,ii,jj)); 2*sqrt(PLKF(4,4,ii,jj))];
    end
end; clear ii jj; % for


%% Plots of Resulting LKF States and Perturbations to Help Debug

% Plot LKF States Results
figure()
sgtitle('LKF State History Typical Simulation')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xLKF(1,:, 1),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xLKF(1,:, 1) + twoSigmasLKF(1,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xLKF(1,:, 1) - twoSigmasLKF(1,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xLKF(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xLKF(2,:, 1) + twoSigmasLKF(2,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xLKF(2,:, 1) - twoSigmasLKF(2,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xLKF(3,:, 1),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xLKF(3,:, 1) + twoSigmasLKF(3,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xLKF(3,:, 1) - twoSigmasLKF(3,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xLKF(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xLKF(4,:, 1) + twoSigmasLKF(4,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xLKF(4,:, 1) - twoSigmasLKF(4,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


%% Plot the LKF State Estimation Errors

% Plot LKF Estimation Errors
figure()
sgtitle('LKF State Estimation Errors Typical Simulation')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xErrorLKF(1,:, 1),'LineWidth',1.5)
ylabel('X Error [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasLKF(1,:, 1), 'r--','LineWidth',1.5)
plot(tvec, -twoSigmasLKF(1,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xErrorLKF(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ Error [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasLKF(2,:, 1), 'r--','LineWidth',1.5)
plot(tvec, -twoSigmasLKF(2,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xErrorLKF(3,:, 1),'LineWidth',1.5)
ylabel('Y Error [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasLKF(3,:, 1), 'r--','LineWidth',1.5)
plot(tvec, -twoSigmasLKF(3,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xErrorLKF(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{Y}$ Error [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasLKF(4,:, 1), 'r--','LineWidth',1.5)
plot(tvec, -twoSigmasLKF(4,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


%% Plot Meas. Innovations to Help Debug

% % Plot Meas Innovations 
% figure()
% sgtitle('Typical Simulation LKF: Measurement Innovations')
% 
% % Rho
% subplot(3,1,1)
% grid on; box on; hold on;
% for jj = 1:nTimesteps
%     perturbInnoVec = eyLKF{jj,1};
%     measInnoVec = delyLKF{jj,1};
%     for ii = 1:3:length(perturbInnoVec)
%         scatter(tvec(jj),perturbInnoVec(ii,1), 'rx')
%         scatter(tvec(jj),measInnoVec(ii,1), 'bx')
%     end
% end; clear ii jj;
% ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
% xlabel('Time [s]')
% legend('Perturbed Innovation', 'Measurement Innovation (deltayk)', 'Location', 'best')
% 
% % Rhodot
% subplot(3,1,2)
% grid on; box on; hold on;
% for jj = 1:nTimesteps
%     perturbInnoVec = eyLKF{jj,1};
%     measInnoVec = delyLKF{jj,1};
%     for ii = 2:3:length(perturbInnoVec)
%         scatter(tvec(jj),perturbInnoVec(ii,1), 'rx')
%         scatter(tvec(jj),measInnoVec(ii,1), 'bx')
%     end
% end; clear ii jj;
% ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
% xlabel('Time [s]')
% legend('Perturbed Innovation', 'Measurement Innovation (deltayk)', 'Location', 'best')
% 
% % phi
% subplot(3,1,3)
% grid on; box on; hold on;
% for jj = 1:nTimesteps
%     perturbInnoVec = eyLKF{jj,1};
%     measInnoVec = delyLKF{jj,1};
%     for ii = 3:3:length(perturbInnoVec)
%         scatter(tvec(jj),perturbInnoVec(ii,1), 'rx')
%         scatter(tvec(jj),measInnoVec(ii,1), 'bx')
%     end
% end; clear ii jj;
% ylabel('$\phi ^i$ [km]', 'Interpreter', 'latex')
% xlabel('Time [s]')
% legend('Perturbed Innovation', 'Measurement Innovation (deltayk)', 'Location', 'best')


%% NEES and NIS Calculations

% Define alpha for these tests
alpha = 0.05; % 5%

% Initialize output epsilons
epsilonx = zeros(nTimesteps, nRuns); % Scalar for each time
epsilony = NaN * ones(nTimesteps, nRuns); % Scalar for each time

% Loop through and calculate these magnitudes using values from the LKF
for jj = 1:nRuns % Loop through each run
    for ii = 1:nTimesteps % Loop through each time

% %         Calculate both magnitudes
%         epsilonx(ii,jj) = exLKF(:,ii,jj)' * (PLKF(:,:,ii,jj))^(-1) * exLKF(:,ii,jj);
%         if ~isempty(eyLKF{ii,jj}) % Make sure there's at least one visible station
%             epsilony(ii,jj) = eyLKF{ii,jj}' * (SLKF{ii,jj})^(-1) * eyLKF{ii,jj};
%         end % if

        % Calculate both magnitudes
        epsilonx(ii,jj) = xErrorLKF(:,ii,jj)' / (PLKF(:,:,ii,jj)) * xErrorLKF(:,ii,jj);
        if ~isempty(yErrorLKF{ii,jj}) % Make sure there's at least one visible station
            epsilony(ii,jj) = yErrorLKF{ii,jj}' / (SLKF{ii,jj}) * yErrorLKF{ii,jj};
        end % if

    end % for ii
end; clear ii jj; % for jj

% Calculate empirical sample averages
epsilonxBar = sum(epsilonx,2) ./ nRuns; % 1401x1
epsilonyBar = sum(epsilony,2) ./ nRuns; % 1401x1

% Determine r1 and r2 bounds for x
r1x = chi2inv(alpha/2, nRuns*length(x0)) ./ nRuns;
r1x = r1x * ones(nTimesteps,1);
r2x = chi2inv(1 - (alpha/2), nRuns*length(x0)) ./ nRuns;
r2x = r2x * ones(nTimesteps,1);


% Determine r1 and r2 bounds for y

% Initialize outputs
r1y = NaN * ones(nTimesteps,1);
r2y = NaN * ones(nTimesteps,1);

% Loop through timesteps
for ii = 1:nTimesteps

    % Determine size of measurement vec at each time
    pTemp = size(ydata{ii},2);

    % Determine if at least one station is visible
    if ~all(isnan(ydata{ii}))
        
        % At least one visible station, calculate r1 and r2
        r1y(ii) = chi2inv(alpha/2, nRuns*pTemp) ./ nRuns;
        r2y(ii) = chi2inv(1 - (alpha/2), nRuns*pTemp) ./ nRuns;

    end % if

end; clear ii; % for


%% Plot NEES Test

% Plot NEES test statistic points 
figure
grid on; box on; hold on;
% Plot the statistics
scatter(tvec(1:nTimesteps), epsilonxBar, 'b')
% Plot the bounds
plot(tvec(1:nTimesteps), r1x, 'r--')
plot(tvec(1:nTimesteps), r2x, 'r--')
title('NEES Test Results: LKF')
xlabel('Time [s]')
ylabel('NEES Statistic, $\bar{\epsilon_x}$', 'Interpreter', 'latex')
legend('NEES @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')
% ylim([min(r1x)-2 max(r2x)+2])


%% Plot NIS Test

% Plot NIS test statistic points 
figure
grid on; box on; hold on;
% Plot the statistics
scatter(tvec(1:nTimesteps), epsilonyBar, 'b')
% Plot the bounds
plot(tvec(1:nTimesteps), r1y, 'r--')
plot(tvec(1:nTimesteps), r2y, 'r--')
title('NIS Test Results: LKF')
xlabel('Time [s]')
ylabel('NIS Statistic, $\bar{\epsilon_y}$', 'Interpreter', 'latex')
legend('NIS @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')
% ylim([min(r1y)-2 max(r2y)+2])


%% Problem 5

%% Tune Initialization of EKF

% After tuning, this P0 was selected 
P0 = diag([1e-10 1e-10 1e-10 1e-10]);

Pplus0 = P0;

% Attempt of initalizing xplus0
xplus0 = x0_pert; % Pert + x0


%% Generate TMT Simulation Data

% Initialize output truth model data
xTMTE = zeros(4, nTimesteps, nRuns);

% Run Monte Carlo Sims and extract TMT data
for ii = 1:nRuns
    xTMTE(:,:,ii) = simulateTruthModel(nTimesteps, Qtrue, x0_pert, P0_noise, deltaT);
end; clear ii; % for


%% Plot Noisy Simulated Ground Truth States for One Instance

% Plot Noisy Simulated Ground Truth States
figure()
sgtitle('Typical Simulation EKF: Noisy Ground Truth States')

% X
subplot(4,1,1)
grid on; box on;
plot(tvec,xTMTE(1,:, 1),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')

% Xdot
subplot(4,1,2)
grid on; box on;
plot(tvec,xTMTE(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Y
subplot(4,1,3)
grid on; box on;
plot(tvec,xTMTE(3,:, 1),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')

% Ydot
subplot(4,1,4)
grid on; box on;
plot(tvec,xTMTE(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')


%% Simulate Noisy Measurements

% Initialize output  data
yNoisyE = NaN * ones(3*nStations, nTimesteps, nRuns);

% Simulate noisy measurements with created function 
for ii = 1:nRuns
    yNoisyE(:,:,ii) = simulateNoisyMeas(yNL(:,1:nTimesteps), Rtrue, nStations);
end; clear ii; % for


%% Plot Noisy Simulated Data

% Plot NL Measurements
figure()
sgtitle('Typical Simulation EKF: Noisy Measurements')

% Rho
subplot(4,1,1)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisyE((ii*3) - 2,:, 1), 'x')
end; clear ii;
ylabel('$\rho ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% Rhodot
subplot(4,1,2)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisyE((ii*3) - 1,:, 1))
end; clear ii;
ylabel('$\dot{\rho} ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% phi
subplot(4,1,3)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec,yNoisyE((ii*3),:, 1), 's')
end; clear ii;
ylabel('$\phi ^i$ [km]', 'Interpreter', 'latex')
xlabel('Time [s]')

% station IDs
subplot(4,1,4)
grid on; box on; hold on;
for ii = 1:nStations
    scatter(tvec, stationIDsMat(ii,1:nTimesteps) , '^')
end; clear ii;
ylabel('Visible Station IDs')
xlabel('Time [s]')


%% Run the EKF

% Initialize outputs
xEKF = zeros(4, nTimesteps, nRuns);
xErrorEKF = zeros(4, nTimesteps, nRuns);
PEKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nRuns); % 4x4x1401xnRuns - Combined so not dependent on stations
SEKF = cell(nTimesteps, nRuns);
yErrorEKF = cell(nTimesteps, nRuns);

% Run EKF
for ii = 1:nRuns
    [xEKF(:,:,ii), xErrorEKF(:,:,ii), PEKF(:,:,:,ii), SEKF(:,ii), yErrorEKF(:,ii)] = runEKFCombined(xTMTE(:,:,ii), yNoisyE(:,:,ii), Qtrue, Rtrue, xplus0, Pplus0, deltaT, RE, omegaE, mu, tvec, nTimesteps, nStations);
end; clear ii; % for


%% Calculate 2sigma Bounds

% Initialize output
twoSigmasEKF = zeros(4, nTimesteps, nRuns);

% Loop through and calculate 2 sigma bounds
for jj = 1:nRuns
    for ii = 1:nTimesteps
        twoSigmasEKF(:, ii, jj) = [2*sqrt(PEKF(1,1,ii,jj)); 2*sqrt(PEKF(2,2,ii,jj)); 2*sqrt(PEKF(3,3,ii,jj)); 2*sqrt(PEKF(4,4,ii,jj))];
    end
end; clear ii jj; % for


%% Plots of Resulting EKF States to Help Debug

% Plot EKF States Results
figure()
sgtitle('EKF State History')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xEKF(1,:, 1),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xEKF(1,:, 1) + twoSigmasEKF(1,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xEKF(1,:, 1) - twoSigmasEKF(1,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xEKF(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xEKF(2,:, 1) + twoSigmasEKF(2,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xEKF(2,:, 1) - twoSigmasEKF(2,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xEKF(3,:, 1),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xEKF(3,:, 1) + twoSigmasEKF(3,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xEKF(3,:, 1) - twoSigmasEKF(3,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xEKF(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xEKF(4,:, 1) + twoSigmasEKF(4,:, 1), 'r--','LineWidth',1.5)
plot(tvec, xEKF(4,:, 1) - twoSigmasEKF(4,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


%% Plot the EKF State Estimation Errors

% Plot EKF Estimation Errors
figure()
sgtitle('EKF State Estimation Errors')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xErrorEKF(1,:, 1),'LineWidth',1.5)
ylabel('X Error [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasEKF(1,:, 1), 'r--','LineWidth',1.5)
plot(tvec, - twoSigmasEKF(1,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')
% ylim([-50 50])

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xErrorEKF(2,:, 1),'LineWidth',1.5)
ylabel('$\dot{X}$ Error [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasEKF(2,:, 1), 'r--','LineWidth',1.5)
plot(tvec, - twoSigmasEKF(2,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')
% ylim([-.5 .5])

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xErrorEKF(3,:, 1),'LineWidth',1.5)
ylabel('Y Error [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasEKF(3,:, 1), 'r--','LineWidth',1.5)
plot(tvec, - twoSigmasEKF(3,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')
% ylim([-50 50])

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xErrorEKF(4,:, 1),'LineWidth',1.5)
ylabel('$\dot{Y}$ Error [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, twoSigmasEKF(4,:, 1), 'r--','LineWidth',1.5)
plot(tvec, - twoSigmasEKF(4,:, 1), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')
% ylim([-.5 .5])


%% NEES and NIS Calculations

% Define alpha for these tests
alphaE = 0.05; % 5%

% Initialize output epsilons
epsilonxE = zeros(nTimesteps, nRuns); % Scalar for each time
epsilonyE = NaN * ones(nTimesteps, nRuns); % Scalar for each time

% Loop through and calculate these magnitudes using values from the LKF
for jj = 1:nRuns % Loop through each run
    for ii = 1:nTimesteps % Loop through each time

        % Calculate both magnitudes
        epsilonxE(ii,jj) = xErrorEKF(:,ii,jj)' / (PEKF(:,:,ii,jj)) * xErrorEKF(:,ii,jj);
        if ~isempty(yErrorEKF{ii,jj}) % Make sure there's at least one visible station
            epsilonyE(ii,jj) = yErrorEKF{ii,jj}' / (SEKF{ii,jj}) * yErrorEKF{ii,jj};
        end % if

    end % for ii
end; clear ii jj; % for jj

% Calculate empirical sample averages
epsilonxBarE = sum(epsilonxE,2) ./ nRuns; % 1401x1
epsilonyBarE = sum(epsilonyE,2) ./ nRuns; % 1401x1

% Determine r1 and r2 bounds for x
r1xE = chi2inv(alpha/2, nRuns*length(x0)) ./ nRuns;
r1xE = r1xE * ones(nTimesteps,1);
r2xE = chi2inv(1 - (alpha/2), nRuns*length(x0)) ./ nRuns;
r2xE = r2xE * ones(nTimesteps,1);


% Determine r1 and r2 bounds for y

% Initialize outputs
r1yE = NaN * ones(nTimesteps,1);
r2yE = NaN * ones(nTimesteps,1);

% Loop through timesteps
for ii = 1:nTimesteps

    % Determine size of measurement vec at each time
    pTemp = size(ydata{ii},2);

    % Determine if at least one station is visible
    if ~all(isnan(ydata{ii}))
        
        % At least one visible station, calculate r1 and r2
        r1yE(ii) = chi2inv(alpha/2, nRuns*pTemp) ./ nRuns;
        r2yE(ii) = chi2inv(1 - (alpha/2), nRuns*pTemp) ./ nRuns;

    end % if

end; clear ii; % for


%% Plot NEES Test

% Plot NEES test statistic points 
figure
grid on; box on; hold on;
% Plot the statistics
scatter(tvec(1:nTimesteps), epsilonxBarE, 'b')
% Plot the bounds
plot(tvec(1:nTimesteps), r1xE, 'r--')
plot(tvec(1:nTimesteps), r2xE, 'r--')
title('NEES Test Results: EKF')
xlabel('Time [s]')
ylabel('NEES Statistic, $\bar{\epsilon_x}$', 'Interpreter', 'latex')
legend('NEES @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')
% ylim([min(r1x)-2 max(r2x)+2])


%% Plot NIS Test

% Plot NIS test statistic points 
figure
grid on; box on; hold on;
% Plot the statistics
scatter(tvec(1:nTimesteps), epsilonyBarE, 'b')
% Plot the bounds
plot(tvec(1:nTimesteps), r1yE, 'r--')
plot(tvec(1:nTimesteps), r2yE, 'r--')
title('NIS Test Results: EKF')
xlabel('Time [s]')
ylabel('NIS Statistic, $\bar{\epsilon_y}$', 'Interpreter', 'latex')
legend('NIS @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')


%% Problem 6

%% Run the LKF

% % Initialize outputs
% xPerturbsLKFData = zeros(4, nTimesteps);
% xLKFData = zeros(4, nTimesteps);
% xErrorLKFData = zeros(4, nTimesteps);
% exLKFData = zeros(4, nTimesteps);
% PLKFData = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps); % 4x4x1401xnRuns - Combined so not dependent on stations
% SLKFData = cell(nTimesteps);
% yErrorLKFData = cell(nTimesteps);
% eyLKFData = cell(nTimesteps);

% Redefine time so that whole orbit occurs
load('orbitdeterm_finalproj_KFdata.mat')
nTimesteps = length(tvec);
% Call ode45 to find nominal trajectory
[tOde,xNom] = ode45(@(t_total,x_total) nonLinearOrbit(t_total,x_total,mu), tspan, x0, tolPts);
% Rearrange output
xNom = interp1(tOde, xNom, tvec);
xNom = xNom';

% Use tuned P0
P0 = diag([1e-10 1e-10 1e-10 1e-10]);

Pplus0 = P0;

% Run LKF
[xPerturbsLKFData, xLKFData, PLKFData, SLKFData] = runLKFwithData(xNom, yNom, ydata, RE, omegaE, tvec, Qtrue, Rtrue, Ftilde, Omegatilde, delxplus0, Pplus0);

%% Calculate 2sigma Bounds

% Initialize output
twoSigmasLKFData = zeros(4, nTimesteps);

% Loop through and calculate 2 sigma bounds
for ii = 1:nTimesteps
    twoSigmasLKFData(:, ii) = [2*sqrt(PLKFData(1,1,ii)); 2*sqrt(PLKFData(2,2,ii)); 2*sqrt(PLKFData(3,3,ii)); 2*sqrt(PLKFData(4,4,ii))];
end; clear ii; % for

%% Plots of Resulting LKF States with Bounds

% Plot LKF States Results
figure()
sgtitle('LKF State History for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xLKFData(1,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xLKFData(1,:) + twoSigmasLKFData(1,:), 'r--','LineWidth',1.5)
plot(tvec,xLKFData(1,:) - twoSigmasLKFData(1,:), 'r--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xLKFData(2,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xLKFData(2,:) + twoSigmasLKFData(2,:), 'r--','LineWidth',1.5)
plot(tvec,xLKFData(2,:) - twoSigmasLKFData(2,:), 'r--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xLKFData(3,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xLKFData(3,:) + twoSigmasLKFData(3,:), 'r--','LineWidth',1.5)
plot(tvec,xLKFData(3,:) - twoSigmasLKFData(3,:), 'r--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xLKFData(4,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xLKFData(4,:) + twoSigmasLKFData(4,:), 'r--','LineWidth',1.5)
plot(tvec,xLKFData(4,:) - twoSigmasLKFData(4,:), 'r--','LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


% Plot LKF States Results Zoomed in
figure()
sgtitle('LKF State History Last 1000s for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec(end-100:end),xLKFData(1,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xLKFData(1,end-100:end) + twoSigmasLKFData(1,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xLKFData(1,end-100:end) - twoSigmasLKFData(1,end-100:end), 'r--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec(end-100:end),xLKFData(2,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xLKFData(2,end-100:end) + twoSigmasLKFData(2,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xLKFData(2,end-100:end) - twoSigmasLKFData(2,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec(end-100:end),xLKFData(3,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xLKFData(3,end-100:end) + twoSigmasLKFData(3,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xLKFData(3,end-100:end) - twoSigmasLKFData(3,end-100:end), 'r--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec(end-100:end),xLKFData(4,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xLKFData(4,end-100:end) + twoSigmasLKFData(4,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xLKFData(4,end-100:end) - twoSigmasLKFData(4,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


%% Run the EKF

% % Initialize outputs
% xEKFData = zeros(4, nTimesteps);
% xErrorEKFData = zeros(4, nTimesteps);
% PEKFData = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps); % 4x4x1401xnRuns - Combined so not dependent on stations
% SEKFData = cell(nTimesteps);
% yErrorEKFData = cell(nTimesteps);

% Pplus0 = diag([1e2, 1e-2, 1e2, 1e-2]);

% Use tuned P0
P0 = diag([1e-10 1e-10 1e-10 1e-10]);

Pplus0 = P0;

% Run EKF
[xEKFData, PEKFData, SEKFData] = runEKFwithData(ydata, Qtrue, Rtrue, xplus0, Pplus0, deltaT, RE, omegaE, mu, tvec, nTimesteps, nStations);


%% Calculate 2sigma Bounds

% Initialize output
twoSigmasEKFData = zeros(4, nTimesteps);

% Loop through and calculate 2 sigma bounds
for ii = 1:nTimesteps
    twoSigmasEKFData(:, ii) = [2*sqrt(PEKFData(1,1,ii)); 2*sqrt(PEKFData(2,2,ii)); 2*sqrt(PEKFData(3,3,ii)); 2*sqrt(PEKFData(4,4,ii))];
end; clear ii; % for


%% Plots of Resulting EKF States with Bounds

% Plot EKF States Results
figure()
sgtitle('EKF State History for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xEKFData(1,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xEKFData(1,:) + twoSigmasEKFData(1,:), 'r--','LineWidth',1.5)
plot(tvec,xEKFData(1,:) - twoSigmasEKFData(1,:), 'r--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xEKFData(2,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xEKFData(2,:) + twoSigmasEKFData(2,:), 'r--','LineWidth',1.5)
plot(tvec,xEKFData(2,:) - twoSigmasEKFData(2,:), 'r--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xEKFData(3,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xEKFData(3,:) + twoSigmasEKFData(3,:), 'r--','LineWidth',1.5)
plot(tvec,xEKFData(3,:) - twoSigmasEKFData(3,:), 'r--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xEKFData(4,:),'LineWidth',1.5)
% Plot bounds
plot(tvec,xEKFData(4,:) + twoSigmasEKFData(4,:), 'r--','LineWidth',1.5)
plot(tvec,xEKFData(4,:) - twoSigmasEKFData(4,:), 'r--','LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


% Plot EKF States Results Zoomed in
figure()
sgtitle('EKF State History Last 1000s for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec(end-100:end),xEKFData(1,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xEKFData(1,end-100:end) + twoSigmasEKFData(1,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xEKFData(1,end-100:end) - twoSigmasEKFData(1,end-100:end), 'r--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec(end-100:end),xEKFData(2,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xEKFData(2,end-100:end) + twoSigmasEKFData(2,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xEKFData(2,end-100:end) - twoSigmasEKFData(2,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec(end-100:end),xEKFData(3,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xEKFData(3,end-100:end) + twoSigmasEKFData(3,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xEKFData(3,end-100:end) - twoSigmasEKFData(3,end-100:end), 'r--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec(end-100:end),xEKFData(4,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xEKFData(4,end-100:end) + twoSigmasEKFData(4,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xEKFData(4,end-100:end) - twoSigmasEKFData(4,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

%% AQ: UKF

%% Tune Initialization of UKF

% After tuning, this P0 was selected
P0 = diag([1e-10 1e-10 1e-10 1e-10]);
Pplus0 = P0;

% Attempt of initalizing xplus0
xplus0 = x0_pert; % Pert + x0


%% Generate TMT Simulation Data

% Define time vector for Monte Carlo sims
nTimesteps = 100; % 100 timesteps corresponds to 1000 s

% Adjust timevector 
tvec = tvec(1:nTimesteps);

% Initialize output truth model data
xTMTU = zeros(4, nTimesteps, nRuns);

% Run Monte Carlo Sims and extract TMT data
for ii = 1:nRuns
    xTMTU(:,:,ii) = simulateTruthModel(nTimesteps, Qtrue, x0_pert, P0_noise, deltaT);
end; clear ii; % for


%% Simulate Noisy Measurements

% Initialize output  data
yNoisyU = NaN * ones(3*nStations, nTimesteps, nRuns);

% Simulate noisy measurements with created function 
for ii = 1:nRuns
    yNoisyU(:,:,ii) = simulateNoisyMeas(yNL(:,1:nTimesteps), Rtrue, nStations);
end; clear ii; % for


%% Run the UKF

% Initialize outputs
xUKF = zeros(4, nTimesteps, nRuns);
xErrorUKF = zeros(4, nTimesteps, nRuns);
PUKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nRuns); % 4x4x1401xnRuns - Combined so not dependent on stations
SUKF = cell(nTimesteps, nRuns);
yErrorUKF = cell(nTimesteps, nRuns);

% Run UKF
for ii = 1:nRuns
    [xUKF(:,:,ii), PUKF(:,:,:,ii), SUKF(:,ii), xErrorUKF(:,:,ii), yErrorUKF(:,ii)] = runUKFCombined(yNoisyU(:,:,ii), Qtrue, Rtrue, xplus0, Pplus0, nTimesteps, deltaT, mu, nStations, tvec, RE, omegaE, xTMTU(:,:,ii));
end; clear ii; % for


%% Calculate 2sigma Bounds

% Initialize output
twoSigmasUKF = zeros(4, nTimesteps, nRuns);

% Loop through and calculate 2 sigma bounds
for jj = 1:nRuns
    for ii = 1:nTimesteps
        twoSigmasUKF(:, ii, jj) = [2*sqrt(PUKF(1,1,ii,jj)); 2*sqrt(PUKF(2,2,ii,jj)); 2*sqrt(PUKF(3,3,ii,jj)); 2*sqrt(PUKF(4,4,ii,jj))];
    end
end; clear ii jj; % for


%% NEES and NIS Calculations

% Define alpha for these tests
alphaU = 0.05; % 5%

% Initialize output epsilons
epsilonxU = zeros(nTimesteps, nRuns); % Scalar for each time
epsilonyU = NaN * ones(nTimesteps, nRuns); % Scalar for each time

% Loop through and calculate these magnitudes using values from the LKF
for jj = 1:nRuns % Loop through each run
    for ii = 1:nTimesteps % Loop through each time

        % Calculate both magnitudes
        epsilonxU(ii,jj) = xErrorUKF(:,ii,jj)' / (PUKF(:,:,ii,jj)) * xErrorUKF(:,ii,jj);
        if ~isempty(yErrorUKF{ii,jj}) % Make sure there's at least one visible station
            epsilonyU(ii,jj) = yErrorUKF{ii,jj}' / (SUKF{ii,jj}) * yErrorUKF{ii,jj};
        end % if

    end % for ii
end; clear ii jj; % for jj

% Calculate empirical sample averages
epsilonxBarU = sum(epsilonxU,2) ./ nRuns; % 1401x1
epsilonyBarU = sum(epsilonyU,2) ./ nRuns; % 1401x1

% Determine r1 and r2 bounds for x
r1xU = chi2inv(alphaU/2, nRuns*length(x0)) ./ nRuns;
r1xU = r1xU * ones(nTimesteps,1);
r2xU = chi2inv(1 - (alphaU/2), nRuns*length(x0)) ./ nRuns;
r2xU = r2xU * ones(nTimesteps,1);


% Determine r1 and r2 bounds for y

% Initialize outputs
r1yU = NaN * ones(nTimesteps,1);
r2yU = NaN * ones(nTimesteps,1);

% Loop through timesteps
for ii = 1:nTimesteps

    % Determine size of measurement vec at each time
    pTemp = size(ydata{ii},2);

    % Determine if at least one station is visible
    if ~all(isnan(ydata{ii}))
        
        % At least one visible station, calculate r1 and r2
        r1yU(ii) = chi2inv(alphaU/2, nRuns*pTemp) ./ nRuns;
        r2yU(ii) = chi2inv(1 - (alphaU/2), nRuns*pTemp) ./ nRuns;

    end % if

end; clear ii; % for

% 
% %% Plot NEES Test
% 
% % Plot NEES test statistic points 
% figure
% grid on; box on; hold on;
% % Plot the statistics
% scatter(tvec(1:nTimesteps), epsilonxBarU, 'b')
% % Plot the bounds
% plot(tvec(1:nTimesteps), r1xU, 'r--')
% plot(tvec(1:nTimesteps), r2xU, 'r--')
% title('NEES Test Results: UKF')
% xlabel('Time [s]')
% ylabel('NEES Statistic, $\bar{\epsilon_x}$', 'Interpreter', 'latex')
% legend('NEES @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')


%% Plot NIS Test

% Plot NIS test statistic points 
figure
grid on; box on; hold on;
% Plot the statistics
scatter(tvec(1:nTimesteps), epsilonyBarU, 'b')
% Plot the bounds
plot(tvec(1:nTimesteps), r1yU, 'r--')
plot(tvec(1:nTimesteps), r2yU, 'r--')
title('NIS Test Results: UKF')
xlabel('Time [s]')
ylabel('NIS Statistic, $\bar{\epsilon_y}$', 'Interpreter', 'latex')
legend('NIS @ Timestep', 'r1 & r2 Bounds', 'Location', 'best')


%% Run the UKF with Data

% Redefine time so that whole orbit occurs
load('orbitdeterm_finalproj_KFdata.mat')
nTimesteps = length(tvec);

% Initialize outputs
% xUKFData = zeros(4, nTimesteps, nRuns);
% PUKFData = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nRuns); % 4x4x1401xnRuns - Combined so not dependent on stations
% SUKFData = cell(nTimesteps, nRuns);

% Run UKF
[xUKFData, PUKFData, SUKFData] = runUKFwithData(ydata, Qtrue, Rtrue, xplus0, Pplus0, nTimesteps, deltaT, mu, nStations, tvec, RE, omegaE);


%% Calculate 2sigma Bounds

% Initialize output
twoSigmasUKFData = zeros(4, nTimesteps);

% Loop through and calculate 2 sigma bounds
for ii = 1:nTimesteps
    twoSigmasUKFData(:, ii) = [2*sqrt(PUKFData(1,1,ii)); 2*sqrt(PUKFData(2,2,ii)); 2*sqrt(PUKFData(3,3,ii)); 2*sqrt(PUKFData(4,4,ii))];
end; clear ii; % for


%% Plots of Resulting UKF States with Data

% Plot UKF States Results
figure()
sgtitle('UKF State History for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec,xUKFData(1,:),'LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xUKFData(1,:) + twoSigmasUKFData(1,:), 'r--','LineWidth',1.5)
plot(tvec, xUKFData(1,:) - twoSigmasUKFData(1,:), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xUKFData(2,:),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xUKFData(2,:) + twoSigmasUKFData(2,:), 'r--','LineWidth',1.5)
plot(tvec, xUKFData(2,:) - twoSigmasUKFData(2,:), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xUKFData(3,:),'LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xUKFData(3,:) + twoSigmasUKFData(3,:), 'r--','LineWidth',1.5)
plot(tvec, xUKFData(3,:) - twoSigmasUKFData(3,:), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xUKFData(4,:),'LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
% Plot bounds
plot(tvec, xUKFData(4,:) + twoSigmasUKFData(4,:), 'r--','LineWidth',1.5)
plot(tvec, xUKFData(4,:) - twoSigmasUKFData(4,:), 'r--','LineWidth',1.5)
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')


% Plot UKF States Results Zoomed in
figure()
sgtitle('UKF State History Last 1000s for Data Log from Canvas')

% X
subplot(4,1,1)
grid on; box on; hold on;
plot(tvec(end-100:end),xUKFData(1,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xUKFData(1,end-100:end) + twoSigmasUKFData(1,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xUKFData(1,end-100:end) - twoSigmasUKFData(1,end-100:end), 'r--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec(end-100:end),xUKFData(2,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xUKFData(2,end-100:end) + twoSigmasUKFData(2,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xUKFData(2,end-100:end) - twoSigmasUKFData(2,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec(end-100:end),xUKFData(3,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xUKFData(3,end-100:end) + twoSigmasUKFData(3,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xUKFData(3,end-100:end) - twoSigmasUKFData(3,end-100:end), 'r--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')

% Ydot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec(end-100:end),xUKFData(4,end-100:end),'LineWidth',1.5)
% Plot bounds
plot(tvec(end-100:end),xUKFData(4,end-100:end) + twoSigmasUKFData(4,end-100:end), 'r--','LineWidth',1.5)
plot(tvec(end-100:end),xUKFData(4,end-100:end) - twoSigmasUKFData(4,end-100:end), 'r--','LineWidth',1.5)
ylabel('$\dot{Y}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('State', '2$\sigma$ Bounds', 'Interpreter', 'latex')



