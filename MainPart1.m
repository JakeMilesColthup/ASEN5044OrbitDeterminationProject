% ASEN 5044
% Billy Weigl, Jake Miles, Adam Buencamino
% Project Part 1

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
perturb = [0; 0.075; 0; -0.021];
x0_pert = x0 + perturb;
% v0 = r0*sqrt(mu/r0^3); % km/s
% omega0 = v0/r0; % rad/s
% ustar = [0; 0]; % Nominal control inputs should be 0 

% % Define nominal station positions - the initial location of only the stations that can see the satellite at the satellite's initial position.
% % All other stations will be empty
% [stationVec, ~] = checkVisibleStations(nStations, RE, omegaE, x0, 0);
% 
% % Linearize system and find the DT linearized model matrices
% [Ftilde0, Gtilde0, Omegatilde0, Htilde0] = computeLinearizedDTModel(x0, deltaT, mu, stationVec);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Might not be necessary
% Note: F, G, and Omega are time invariant, but H is time-varying.

% % Calculate stability of the time invariant system
% eigValuesLTI = eig(Ftilde);
% 
% % Calculate the observability matrix
% obsMat = obsv(Ftilde,Htilde0);
% 
% % Determine observability
% observable = rank(obsMat);
% if observable == size(obsMat,2)
%     fprintf('LTI System is Observable! \n')
% end

% Make time vector for just one orbital period
% onePeriod = 2*pi*sqrt(r0^3/mu); % s
% tvecOnePeriod = 0:10:onePeriod; % s


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

% Simulate linearized DT dynamics and measurement models
for ii = 1:length(tvec)-1

    % Linearize system and find F, G, Omega
    [Ftilde(:,:,ii+1), Gtildek, Omegatilde(:,:,ii+1)] = computeLinearizedDyn(xNom(:,ii), deltaT, mu); 

    % Simulate DT Dynamics
    xLinearized(:,ii+1) = xNom(:,ii+1) + Ftilde(:,:,ii+1)*( xLinearized(:,ii) - xNom(:,ii) ) + Gtildek*[0;0] + Omegatilde(:,:,ii+1)*[0;0]; % Note: no control input perturbations, no process noise

    % Calculate information of the visible stations
    [stationVec, stationIDs] = checkVisibleStations(nStations, RE, omegaE, xNL(:,ii+1), tvec(ii+1));

    % Calculate H
    Htilde(:,:,ii+1) = computeLinearizedH(xNom(:,ii+1), stationVec);

    % Simulate measurement model
    yNL(:,ii+1) = computeYNL(nStations, stationIDs, RE, omegaE, xNL(:,ii+1), tvec(ii+1));
    yNom(:,ii+1) = computeYNL(nStations, stationIDs, RE, omegaE, xNom(:,ii+1), tvec(ii+1));
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

% Xdot
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
plot(tvec,xLinearized(1,:),'--','LineWidth',1.5)
ylabel('X [km]')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Xdot
subplot(4,1,2)
grid on; box on; hold on;
plot(tvec,xNL(2,:),'LineWidth',1.5)
plot(tvec,xLinearized(2,:),'--','LineWidth',1.5)
ylabel('$\dot{X}$ [km/s]', 'Interpreter', 'latex')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Y
subplot(4,1,3)
grid on; box on; hold on;
plot(tvec,xNL(3,:),'LineWidth',1.5)
plot(tvec,xLinearized(3,:),'--','LineWidth',1.5)
ylabel('Y [km]')
xlabel('Time [s]')
legend('Full Nonlinear', 'Linearized', 'Location', 'best')

% Xdot
subplot(4,1,4)
grid on; box on; hold on;
plot(tvec,xNL(4,:),'LineWidth',1.5)
plot(tvec,xLinearized(4,:),'--','LineWidth',1.5)
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

% Xdot
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


