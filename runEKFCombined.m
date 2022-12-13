function [xEKF, xErrorEKF, PEKF, SEKF, yErrorEKF] = runEKFCombined(xTMT, yNoisy, Qtrue, Rtrue, xplus0, Pplus0, deltaT, RE, omegaE, mu, tvec, nTimesteps, nStations)
% Run EKF once for all timesteps. Simulates one complete simulation based on simulated noisy measurements.
% Inputs are nominal trajectory (NOT used in the EKF itself, just to calculate error after), noisy measurements, Qtrue and Rtrue from Canvas data,
% nominal state history, initial state and cov of EKF, Earth parameters, time vector from data file, nTimesteps, nStations
% Outputs are state perturbations from the nominal trajectory, total states predicted by EKF, perturbations in measurements from nominal measurements,
% estimated state errors from EKF. Note: xError is ex, yError is ey. PEKF and SEKF are the P and S used for NEES/NIS

%% Tune Noise 
% Tune Q
Q = Qtrue;

% After tuning, this Q was selected
Q = 1e-12 .*eye(2);

% Set R (also assuming R stays constant)
R = Rtrue;
% R = R * 2.5; % for 10 runs
R = R * 2.7; % for 25

%% Initialize Variables

% Determine number of states
n = size(xplus0,1);

% Initialize state estimate outputs
xEKF = zeros(size(xplus0,1), nTimesteps); 

xEKF(:,1) = xplus0;

% Initialize storing of the combined covariance matrices (combined for all stations)
PEKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);
PEKF(:,:,1) = Pplus0;
PEKFMinus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);

% Initialize Sk & ey (Changes size based on how many visible stations there are)
SEKF = cell(nTimesteps,1);
SEKF{1} = R;
yErrorEKF = cell(nTimesteps,1);
yErrorEKF{1} = [];
etildey = cell(nTimesteps,1);
etildey{1} = [];

% Initialize jacobians
Ftilde = zeros(n, n, nTimesteps);
Omegatilde = zeros(n, 2, nTimesteps);
Htilde = NaN * ones(3*nStations, n, nTimesteps);

[Ftilde(:,:,1), ~, Omegatilde(:,:,1)] = computeLinearizedDyn(xEKF(:,1), deltaT, mu); 
[stationVec, ~] = checkVisibleStations(nStations, RE, omegaE, xEKF(:,1), tvec(1));
Htilde(:,:,1) = computeLinearizedH(xEKF(:,1), stationVec);

% Define ode tolerance
tolPts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% EKF Loop

% Loop through all timesteps and implement the EKF to estimate perturbation states
for ii = 2:nTimesteps
 
    % Increase process noise Qkf at 700 seconds into the sim
    if ii == 71
        Q = Q * 10;
    end % if

    % TIME UPDATE STEP FOR TIME K+1 (except Pkp1Minus, that's done later)
    
    % Note: here, k+1 is technically equal to ii-1 because of matlab 1 indexing, but for this purpose it doesn't matter since we will be indexing
    % into vectors with "kp1" which is the element corresponding to the proper k+1. For example, on the first timestep k+1 is ACTUALLY 1, but index
    % 2 into the vectors represents k=1 since Matlab starts indexing at 1 but the first element is for k=0. Therefore, really k = ii-2.
    kp1 = ii;
    k = kp1 - 1;

    % Compute delxkp1Minus using the nonlinear dynamics (approximated with ode45)
    [~,xhatkp1Minus_Dyn] = ode45(@(t_total,x_total) nonLinearOrbit(t_total,x_total,mu), [0 deltaT], xEKF(:,k), tolPts);

    % Store (note: comes out as row vector)
    xhatkp1Minus = xhatkp1Minus_Dyn(end,:)';

    % Note: using nominal control input values so u, ustar, deltau, are all 0 so those parts of the equations go away

    % Find F, Omega
    [Ftilde(:,:,kp1), ~, Omegatilde(:,:,kp1)] = computeLinearizedDyn(xEKF(:,k), deltaT, mu); 

    % MEASUREMENT UPDATE STEP FOR TIME K+1

    % Initialize xhatkp1Plus to current value of xhatkp1Minus
    xhatkp1Plus = xhatkp1Minus;

    % Calculate information of the visible stations
%     [stationVec, visibleStationIds] = checkVisibleStations(nStations, RE, omegaE, xhatkp1Minus, tvec(kp1));
    [stationVec, visibleStationIds] = checkVisibleStationsNoisyMeas(nStations, RE, omegaE, tvec(kp1), yNoisy(:,kp1));

    % Calculate H
    Htilde(:,:,kp1) = computeLinearizedH(xhatkp1Minus, stationVec);

    % Compute yhatkp1Minus
    yhatkp1Minus = computeYNL(nStations, visibleStationIds, RE, omegaE, xhatkp1Minus, tvec(kp1));

    % Initialize stationTracker
    stationTracker = 1;

    % Initialize visibleStationjjs
    visibleStationIds = [];
    visibleStationjjs = [];

    % Loop through each station
    for jj = 1:3:3*nStations

        % Check if this station is visible
        if all(~isnan(yNoisy(jj:jj+2, kp1)))

            % Store this station as being visible
            visibleStationIds = [visibleStationIds, stationTracker];
            visibleStationjjs = [visibleStationjjs, jj];

        end % if

        % Iterate station tracker
        stationTracker = stationTracker + 1;

    end % for

    % Determine if any measurement data is present
    if ~isempty(visibleStationIds)

        % Measurements are present, continue with meas. update

        % Calculate next state and Calculate statstics for tests 

        % Calculate the combined Pkp1Minus
        PEKFMinus(:,:,kp1) = Ftilde(:,:, k) * PEKF(:,:,k) * Ftilde(:,:, k)' + Omegatilde(:,:, k) * Q * Omegatilde(:,:, k)';

        % Loop through visible stations, combine their K and Hs
        Hcombined = [];
        for kk = 1:length(visibleStationIds)
            Hcombined = [Hcombined; Htilde(visibleStationjjs(kk):visibleStationjjs(kk)+2, :, kp1)];
        end

        % Create block diagonal big R matrix for this number of visible stations
        Rbig = [];
        for qq = 1:length(visibleStationjjs)
            Rbig = blkdiag(Rbig, R);
        end

        % Calculate Kcombined
        Kcombined = PEKFMinus(:,:,kp1) * Hcombined' / ( Hcombined * PEKFMinus(:,:,kp1) * Hcombined' + Rbig );

        % Calculate the combined Pkp1Plus
        PEKF(:,:,kp1) = (eye(n) -  Kcombined * Hcombined) * PEKFMinus(:,:,kp1);

        % Calculate SEKF
        SEKF{kp1} = ( Hcombined * PEKFMinus(:,:,kp1) * Hcombined' ) + Rbig;

        % Extract y data into one stacked vec at this timestep
        [yHatVec, yNoisyVec] = reformatYNomWithNoisy(yhatkp1Minus, yNoisy(:,kp1));

        % Ensure that size of yNoisy here matches size of y combined
        assert(length(yNoisyVec) == size(Hcombined,1), 'ERROR: Size of yNoisy does not match size of Hcombined')

        % Calculate NL meas. innovation (actual data minus predicted)
        etildeykp1 = yNoisyVec - yHatVec;

        % Calculate delxkp1Plus with measurement data
        xhatkp1Plus = xhatkp1Plus + Kcombined * etildeykp1;

        % Store meas. innovation
        etildey{kp1} = [etildey{kp1}; etildeykp1];

        % Calculate yError (or usual ey,k)
        yErrorEKF{kp1} = etildey{kp1};


    else % No visible stations, don't measurement update

        % Set P
        PEKF(:,:,kp1) = PEKF(:,:,k);

        % Set Sk = R
        SEKF{kp1} = R;

        % Set eyk to empty
        yErrorEKF{kp1} = [];

        % Don't update any measurement to new xEKF value
        xEKF(:,kp1) = xhatkp1Minus;

    end % if

    % Store calculated state
    xEKF(:,kp1) = xhatkp1Plus;
    
end; clear ii jj kk qq; % for

% Calculate estimated state error using the ground truth data
xErrorEKF = xTMT - xEKF; % AKA usual ex

end % function