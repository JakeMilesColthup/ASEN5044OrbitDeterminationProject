function [xUKF, PUKF, SUKF, xErrorUKF, yErrorUKF] = runUKFCombined(yNoisy, Qtrue, Rtrue, xplus0, Pplus0, nTimesteps, deltaT, mu, nStations, tvec, RE, omegaE, xTMT)
% Run UKF once for all timesteps. Simulates one complete simulation based on measurements.
% Inputs are noisy meas, Qtrue, Rtrue, initial state and cov, nTimesteps, ..., truth model data
% Outputs are state perturbations from the nominal trajectory, total states predicted by UKF, perturbations in measurements from nominal measurements,
% estimated state errors from UKF. Note: xError is ex, yError is ey. PUKF and SUKF are the P and S used for NEES/NIS

%% Tune Noise 
% Tune Q
Q = Qtrue;

% After tuning, this Qkf was selected
Q = 1e-15 .*eye(2);

Qbig = diag([0, Q(1,1), 0, Q(2,2)]);

% Note: assuming that Q stays constant throughout sim so Qk = Q

% Set R (also assuming R stays constant)
R = Rtrue;
R = R / 2.5; % for 25

%% Initialize Variables

% Determine number of states
n = size(xplus0,1);

% Initialize state estimate outputs
xUKF = zeros(size(xplus0,1), nTimesteps); 

xUKF(:,1) = xplus0;

% Initialize storing of the combined covariance matrices (combined for all stations)
PUKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);
PUKF(:,:,1) = Pplus0;
PUKFMinus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);

% Initialize Sk & ey (Changes size based on how many visible stations there are)
SUKF = cell(nTimesteps,1);
SUKF{1} = R;
yErrorUKF = cell(nTimesteps,1);
yErrorUKF{1} = [];

% Initialize storing of other cov matrices
Pxy = cell(nTimesteps,1);
Pxy{1} = [];
Pyy = cell(nTimesteps,1);
Pyy{1} = [];

% Define UKF parameters
kappa = 0;
beta = 2;
alpha = 1.5;
lambda = alpha^2 * (n+kappa) - n;

% Initialize wm and wc
wm = zeros(2*n+1,1);
wc = zeros(2*n+1,1);

% Calculate both
wm(1) = lambda/(n+lambda);
wc(1) = lambda/(n+lambda) + 1 - alpha^2 + beta;
for zz = 2:2*n+1
    wm(zz) = 1/(2*(n+lambda));
    wc(zz) = wm(zz);
end; clear zz; % for

% Define ode tolerance
tolPts = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% UKF Loop

% Loop through all timesteps and implement the UKF to estimate perturbation states
for ii = 2:nTimesteps

    % TIME UPDATE STEP FOR TIME K+1 (except Pkp1Minus, that's done later)
    
    % Note: here, k+1 is technically equal to ii-1 because of matlab 1 indexing, but for this purpose it doesn't matter since we will be indexing
    % into vectors with "kp1" which is the element corresponding to the proper k+1. For example, on the first timestep k+1 is ACTUALLY 1, but index
    % 2 into the vectors represents k=1 since Matlab starts indexing at 1 but the first element is for k=0. Therefore, really k = ii-2.
    kp1 = ii;
    k = kp1 - 1;

    % Find Sk
    Sk = chol(PUKF(:,:,k));

    % Generate sigma pts (9 of them = 2n+1)
    sigmaPtsk = zeros(n,2*n+1); % Initialization
    % Generate the first one
    sigmaPtsk(:,1) = xUKF(:,k);
    % Generate the rest
    for zz = 1+1:n+1
        sigmaPtsk(:,zz) = xUKF(:,k) + (sqrt(n+lambda)) * Sk(zz-1,:)';
    end; clear zz; % for
    for zz = n+1+1:2*n+1
        sigmaPtsk(:,zz) = xUKF(:,k) + (sqrt(n+lambda)) * Sk(zz-n-1,:)';
    end; clear zz; % for

    % Initialize resultant pts
    resultantPts = zeros(n,2*n+1);

    % Compute resultant points using the nonlinear dynamics (approximated with ode45)
    for zz = 1:2*n+1
        [~,resultantDyn] = ode45(@(t_total,x_total) nonLinearOrbit(t_total,x_total,mu), [0 deltaT], sigmaPtsk(:,zz), tolPts);
        resultantPts(:,zz) = resultantDyn(end,:)';
    end; clear zz; % for
    
    % Combine resultant pts to get mean
    sum = zeros(n,1);
    for zz = 1:2*n+1
         sum = sum + wm(zz)*resultantPts(:,zz);
    end; clear zz; % for
    xhatkp1Minus = sum;

    % Combine resultant pts to get cov
    sum = zeros(n,n);
    for zz = 1:2*n+1
         sum = sum + wc(zz) * ( resultantPts(:,zz)-xhatkp1Minus ) * ( resultantPts(:,zz)-xhatkp1Minus )' + Qbig;
    end; clear zz; % for
    PUKFMinus(:,:,kp1) = sum;


    % MEASUREMENT UPDATE STEP FOR TIME K+1

    % Initialize xhatkp1Plus to current value of xhatkp1Minus
    xhatkp1Plus = xhatkp1Minus;

    % Initialize stationTracker
    stationTracker = 1;

    % Initialize visibleStationIds
    visibleStationIds = [];
    visibleStationjjs = [];

    % Loop through each station
    for jj = 1:3:3*nStations

        % Check if this station is visible
        if all(~isnan(yNoisy(jj:jj+2, kp1)))

            % It is visible, store ids

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

        % Find Skp1bar
        Skp1bar = chol(PUKFMinus(:,:,kp1));

        % Generate sigma pts (9 of them = 2n+1)
        sigmaPtskp1 = zeros(n,2*n+1); % Initialization
        % Generate the first one
        sigmaPtskp1(:,1) = xhatkp1Minus;
        % Generate the rest
        for zz = 1+1:n+1
            sigmaPtskp1(:,zz) = xhatkp1Minus + (sqrt(n+lambda)) * Skp1bar(zz-1,:)';
        end; clear zz; % for
        for zz = n+1+1:2*n+1
            sigmaPtskp1(:,zz) = xhatkp1Minus + (sqrt(n+lambda)) * Skp1bar(zz-n-1,:)';
        end; clear zz; % for

        % Extract y measurements into one stacked vec at this timestep
        yNoisyVec = reformatYNom(yNoisy(:,kp1));

        % Initialize resultant pts gamma
        gammaPts = zeros(length(yNoisyVec),2*n+1);

        % Extract y measurements into 36x1 vec with nans
        yNoisyVec36 = yNoisy(:,kp1);

        % Calculate information of the visible stations
        [~, visibleStationIds] = checkVisibleStationsNoisyMeas(nStations, RE, omegaE, tvec(kp1), yNoisyVec36);

        % Compute resultant points using the nonlinear dynamics (approximated with ode45)
        for zz = 1:2*n+1
            gammaPts36 = computeYNL(nStations, visibleStationIds, RE, omegaE, sigmaPtskp1(:,zz), tvec(kp1));
            gammaPts(:,zz) = reformatYNom(gammaPts36);
        end; clear zz; % for

        % Make sure sizes are the same
        assert(length(yNoisyVec) == size(gammaPts,1), 'ERROR: Size of yData does not match size of gammaPts')

        % Combine resultant pts to get mean
        sum = zeros(length(yNoisyVec),1);
        for zz = 1:2*n+1
            sum = sum + wm(zz)*gammaPts(:,zz);
        end; clear zz; % for
        yhatkp1Minus = sum;

        % Create block diagonal big R matrix for this number of visible stations
        Rbig = [];
        for qq = 1:length(visibleStationIds)
            Rbig = blkdiag(Rbig, R);
        end

        % Combine resultant pts to get cov
        sum = zeros(length(yNoisyVec),length(yNoisyVec));
        for zz = 1:2*n+1
            sum = sum + wc(zz) * ( gammaPts(:,zz)-yhatkp1Minus ) * ( gammaPts(:,zz)-yhatkp1Minus )' + Rbig;
        end; clear zz; % for
        Pyy{kp1} = sum;

        % Get Pxy (nxp)
        sum = zeros(n,length(yNoisyVec));
        for zz = 1:2*n+1
            sum = sum + wc(zz) * ( resultantPts(:,zz)-xhatkp1Minus ) * ( gammaPts(:,zz)-yhatkp1Minus )';
        end; clear zz; % for
        Pxy{kp1} = sum;

        % Get Kalman gain matrix (nxp)
        Kkp1 = Pxy{kp1} / (Pyy{kp1});

        % Extract y data into one stacked vec at this timestep
        yHatVec = yhatkp1Minus;

        % Calculate new state
        xhatkp1Plus = xhatkp1Minus + Kkp1 * (yNoisyVec - yHatVec);

        % Calculate Pkp1plus
        PUKF(:,:,kp1) = PUKFMinus(:,:,kp1) - Kkp1 * Pyy{kp1} * Kkp1';

        % Calculate yError (or usual ey,k)
        yErrorUKF{kp1} = yNoisyVec - yHatVec;

        % Calculate S
        SUKF{kp1} = Pyy{kp1};


    else % No visible stations, don't measurement update

        % Set P
        PUKF(:,:,kp1) = PUKF(:,:,k);

        % Set Sk = R
        SUKF{kp1} = R;

        % Set eyk to empty
        yErrorUKF{kp1} = [];

        % Don't update any measurement to new xUKF value
        xUKF(:,kp1) = xhatkp1Minus;

    end % if

    % Store calculated state
    xUKF(:,kp1) = xhatkp1Plus;
    
end; clear ii jj kk qq; % for

% Calculate estimated state error using the ground truth data
xErrorUKF = xTMT - xUKF; % AKA usual ex

end % function