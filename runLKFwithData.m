function [xPerturbsLKF, xLKF, PLKF, SLKF] = runLKFwithData(xNom, yNom, ydata, RE, omegaE, tvec, Qtrue, Rtrue, Ftilde, Omegatilde, delxplus0, Pplus0)
% Run LKF once for all timesteps. Simulates one complete simulation based on nominal traj and measurement data from Canvas.
% Inputs are nominal trajectory data, nominal measurement data, ydata, Qtrue and Rtrue from Canvas data,
% nominal state history, Ftilde for each timestep, Omegatilde for each timestep, Htilde for each timstep, initial state and cov of LKF
% Outputs are state perturbations from the nominal trajectory, total states predicted by LKF, perturbations in measurements from nominal measurements,
% estimated state errors from LKF. Note: xError is ex, yError is ey. PLKF and SLKF are the P and S used for NEES/NIS
% Update: added new ex and ey variables to calculate perturbations from nominal traj not total states.

%% Tune Noise
% Tune Q
Q = Qtrue;
% Q = Q .* 1e7;
% Q = eye(2); % this caused spike starting at ii = 165
% Q = 1e-20 .*eye(2);
% Q = [1e-2 0; 0 1];

% After tuning, this was selected
Q = 1e-12 .*eye(2);

% Set R (assuming R stays constant)
R = Rtrue;
% R = R * 2.5; % for 10 runs
R = R * 2.7; % for 25

%% Initialize Variables

% Determine number of timesteps
nTimesteps = size(xNom,2);

% Determine number of states
n = size(xNom,1);

% Determine number of stations
nStations = size(yNom,1)/3;

% Initialize state estimate outputs
xPerturbsLKF = zeros(size(xNom,1), nTimesteps); 

xPerturbsLKF(:,1) = delxplus0; 

% Initialize storing of the combined covariance matrices (combined for all stations)
PLKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);
PLKF(:,:,1) = Pplus0;
PLKFMinus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);

% Initialize Sk & ey (Changes size based on how many visible stations there are)
SLKF = cell(nTimesteps,1);
SLKF{1} = R;
% yErrorLKF = cell(nTimesteps,1);
% yErrorLKF{1} = [];
% eyLKF = cell(nTimesteps,1);
% eyLKF{1} = [];
% delyLKF = cell(nTimesteps,1);
% delyLKF{1} = [];

% Initialize jacobian
Htilde = NaN * ones(3*nStations, n, nTimesteps);

yDataVec36 = reformatYDataCell36(ydata{1}, nStations);
[stationVec, ~] = checkVisibleStationsNoisyMeas(nStations, RE, omegaE, tvec(1), yDataVec36);
Htilde(:,:,1) = computeLinearizedH(xNom(:,1), stationVec);


%% LKF Loop

% Loop through all timesteps and implement the LKF to estimate perturbation states
for ii = 2:nTimesteps

    % TIME UPDATE STEP FOR TIME K+1 (except Pkp1Minus, that's done later)
    
    % Note: here, k+1 is technically equal to ii-1 because of matlab 1 indexing, but for this purpose it doesn't matter since we will be indexing
    % into vectors with "kp1" which is the element corresponding to the proper k+1. For example, on the first timestep k+1 is ACTUALLY 1, but index
    % 2 into the vectors represents k=1 since Matlab starts indexing at 1 but the first element is for k=0. Therefore, really k = ii-2.
    kp1 = ii;
    k = kp1 - 1;

    % Compute delxkp1Minus
    delxkp1Minus = Ftilde(:,:, k) * xPerturbsLKF(:,k); 

    % Note: using nominal control input values so u, ustar, deltau, are all 0 so those parts of the equations go away


    % MEASUREMENT UPDATE STEP FOR TIME K+1

    % Initialize delxkp1Plus to current value of delxkp1Minus
    delxkp1Plus = delxkp1Minus;

    % Adjust yNom to match ydata
    % Extract y data into one stacked vec at this timestep
    yDataVec = reformatYDataCell(ydata{kp1});

    % Extract y data into 36x1 vec with nans
    yDataVec36 = reformatYDataCell36(ydata{kp1}, nStations);

    % Extract yNom
    [yNomVec, ~] = reformatYNomWithNoisy(yNom(:,kp1), yDataVec36);

    % Recalculate H with this nominal data
    [stationVec, ~] = checkVisibleStationsNoisyMeas(nStations, RE, omegaE, tvec(kp1), yDataVec36);
    Htilde(:,:,kp1) = computeLinearizedH(xNom(:,kp1), stationVec);

    % Loop through each visible station, and then do whole measurement step that many times

    % Initialize stationTracker
    stationTracker = 1;

    % Initialize visibleStationIds
    visibleStationIds = [];
    visibleStationjjs = [];

    % Loop through each station
    for jj = 1:3:3*nStations    

        % Check if this station is visible
        if (~isempty(ydata{kp1}))
            if any(ydata{kp1}(4,:) == stationTracker)

                % It is visible, store ids

                % Store this station as being visible
                visibleStationIds = [visibleStationIds, stationTracker];
                visibleStationjjs = [visibleStationjjs, jj];

            end % if
        end % if

        % Iterate station tracker
        stationTracker = stationTracker + 1;

    end % for

    % Calculate next perturbation and Store current PLKF and SLKF
    if ~isempty(visibleStationIds)
        % Calculate the combined Pkp1Minus
        PLKFMinus(:,:,kp1) = Ftilde(:,:, k) * PLKF(:,:,k) * Ftilde(:,:, k)' + Omegatilde(:,:, k) * Q * Omegatilde(:,:, k)';

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
        Kcombined = PLKFMinus(:,:,kp1) * Hcombined' / ( Hcombined * PLKFMinus(:,:,kp1) * Hcombined' + Rbig );

        % Calculate the combined Pkp1Plus
        PLKF(:,:,kp1) = (eye(n) -  Kcombined * Hcombined) * PLKFMinus(:,:,kp1);

        % Calculate SLKF
        SLKF{kp1} = ( Hcombined * PLKFMinus(:,:,kp1) * Hcombined' ) + Rbig;

        % Ensure that size of ydata here matches size of y combined
        assert(length(yDataVec) == size(Hcombined,1) && all(all(~isnan(Hcombined))), 'ERROR: Size of ydata does not match size of Hcombined')

        % Calculate delykp1
        delykp1 = yDataVec - yNomVec;

        % Calculate yError (or delhatykp1)
%         yErrorLKF{kp1} = yDataVec - ( (Hcombined * delxkp1Minus) + yNomVec );
        delhatykp1 = (Hcombined * delxkp1Minus);

        % Calcualte delxp1plus
        delxkp1Plus = delxkp1Minus + Kcombined * (delykp1 - delhatykp1);

        % Calculate new perturbation ey,k
%         eyLKF{kp1} = (yDataVec - yNomVec) - (Hcombined * delxkp1Minus);

        % Store meas. innovation
%         delyLKF{kp1} = delykp1;

    else % No visible stations

        % Set P
        PLKF(:,:,kp1) = PLKF(:,:,k);

        % Set Sk = R
        SLKF{kp1} = R;

        % Set eyk to empty
%         yErrorLKF{kp1} = [];
%         eyLKF{kp1} = [];
%         delyLKF{kp1} = [];

    end % if

    % Store perturbation calculated
    xPerturbsLKF(:,kp1) = delxkp1Plus; 

end; clear ii jj kk qq; % for

% Calculate total states
xLKF = xNom + xPerturbsLKF;

% Calculate estimated state error using the ground truth data
% xErrorLKF = xTMT - xLKF; % AKA usual ex
% exLKF = (xTMT - xNom) - (xLKF - xNom); % Perturbation ex
% xErrorLKF = xTMT - xLKF;
% exLKF = -1 * (xLKF - xTMT);

end % function