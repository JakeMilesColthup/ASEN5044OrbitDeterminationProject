function xTruth = simulateTruthModelWholeSim(nTimesteps, Qtrue, mu0, P0, tvec)
% Simulates truth model by adding process noise to a randomly perturbated trajectory.
% Input parameters nTimesteps and alpha are specified in the Main script.
% Input Qtrue given in provided data file.
% Input mu0 = nominal x0 from project description.

%% Initialize Sample Perturbated Trajectory

% Sample p(x(0)) to create new initial condition
x0_sample = mvnrnd(mu0, P0);

% Convert to proper 4x1 column vector
x0_sample = x0_sample';

% Initialize output
xTruth = zeros(4, nTimesteps);
xTruth(:,1) = x0_sample;

%% Parameter Definition

% Define parameters
mu = 398600; % km^3/s^2

% Define ode tolerance
tolPts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% Convert Q Provided to Sigma Used for Sampling

% Note: this works because Qtrue is a diagonal matrix that only affects acceleration
sigMat = diag([0, Qtrue(1,1), 0, Qtrue(2,2)]);

%% Truth Model Simulation

% Generate process noise for this timestep
wkMat = mvnrnd(zeros(length(tvec), length(x0_sample)), sigMat);

% Convert to column vector
wkMat = wkMat';

% Call ode45 to push through nonlinear dynamic model
[~,xk_Dyn] = ode45(@(t_total,x_total) nonLinearOrbitNoiseWholeSim(t_total,x_total,mu, wkMat, tvec), tvec, xTruth(:,1), tolPts);

% Add noise to output and store (note: comes out as row vectors)
xTruth(:,2:end) = xk_Dyn(2:end,:)';


end % function