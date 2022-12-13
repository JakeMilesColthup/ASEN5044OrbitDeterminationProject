function [xPerturbsLKF, xLKF, xErrorLKF, PLKF, SLKF, yErrorLKF, exLKF, eyLKF] = runLKF(xTMT, xNom, yNom, yNoisy, Qtrue, Rtrue, Ftilde, Omegatilde, Htilde, delxplus0, Pplus0)
% Run LKF once for all timesteps. Simulates one complete simulation based on nominal traj and simulated measurment data.
% Inputs are truth model data, nominal trajectory data, nominal measurement data, simulated noisy measurement data, Qtrue and Rtrue from Canvas data,
% nominal state history, Ftilde for each timestep, Omegatilde for each timestep, Htilde for each timstep, initial state and cov of LKF
% Outputs are state perturbations from the nominal trajectory, total states predicted by LKF, perturbations in measurements from nominal measurements,
% estimated state errors from LKF. Note: xError is ex, yError is ey. PLKF and SLKF are the P and S used for NEES/NIS
% Note: Chose to sequentially perform the measurement update for each visible station, but can also combine it and do it all at once.
% Update: added new ex and ey variables to calculate perturbations from nominal traj not total states.

%% Tune Noise (EDIT THIS!!!!!!!!)
% Tune Q
Q = Qtrue;
% Q = Q .* 1e4;

% Note: assuming that Q stays constant throughout sim so Qk = Q

% Set R (also assuming R stays constant)
R = Rtrue;

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

% Initialize previous iteration variables to the initial condition
% Note: Chose to make Pkp1Plus a 4D matrix so that each station can keep track of its own P matrix at every timestep
Pkp1Plus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nStations);
for ii = 1:nStations
    Pkp1Plus(:,:,1,ii) = Pplus0;
end; clear ii;
Pkp1Minus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps, nStations);

% Initialize storing of the combined covariance matrices (combined for all stations)
Kkp1Mat = zeros(4,3, nTimesteps, nStations);
PLKF = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);
PLKF(:,:,1) = Pplus0;
PLKFMinus = zeros(size(Pplus0,1), size(Pplus0,2), nTimesteps);

% Initialize Sk & ey (Changes size based on how many visible stations there are)
SLKF = cell(nTimesteps,1);
SLKF{1} = R;
yErrorLKF = cell(nTimesteps,1);
yErrorLKF{1} = [];
eyLKF = cell(nTimesteps,1);
eyLKF{1} = [];

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

    % Extract ydata
    yNoisyVec = yNoisy(:,kp1);

    % Compute delykp1
    delykp1 = yNoisyVec - yNom(:,kp1);

    % Loop through each visible station, and then do whole measurement step that many times

    % Initialize stationTracker
    stationTracker = 1;

    % Initialize visibleStationIds
    visibleStationIds = [];
    visibleStationjjs = [];

    % Loop through each station
    for jj = 1:3:3*nStations    

        % Compute Pkp1Minus
        Pkp1Minus(:,:,kp1,stationTracker) = Ftilde(:,:, k) * Pkp1Plus(:,:,k,stationTracker) * Ftilde(:,:, k)' + Omegatilde(:,:, k) * Q * Omegatilde(:,:, k)';

        % Extract this station's information from Htilde
        Htildei = Htilde(jj:jj+2, :, kp1);

        % Compute Kkp1
        Kkp1 = Pkp1Minus(:,:,kp1,stationTracker) * Htildei' * ( Htildei*Pkp1Minus(:,:,kp1,stationTracker)*Htildei' + R )^(-1);
%         Kkp1 = Pkp1Minus(:,:,kp1,stationTracker) * Htildei' / ( Htildei*Pkp1Minus(:,:,kp1,stationTracker)*Htildei' + R );

        % Store this value of Kkp1
        Kkp1Mat(:,:,kp1,stationTracker) = Kkp1;

        % Check if this station is visible
        if all(~isnan(yNom(jj:jj+2, kp1)))

            % It is visible, perform measurement update

            % Compute Pkp1Plus for next iteration
            Pkp1Plus(:,:,kp1,stationTracker) = (eye(n) -  Kkp1 * Htildei) * Pkp1Minus(:,:,kp1,stationTracker);

            % Calculate delxkp1Plus with measurement data
            delxkp1Plus = delxkp1Plus + Kkp1 * ( delykp1(jj:jj+2) - Htildei * delxkp1Plus );

            % Store this station as being visible
            visibleStationIds = [visibleStationIds, stationTracker];
            visibleStationjjs = [visibleStationjjs, jj];

        else
            
            % It is not visible, prevent Pkp1Plus from being stored as NaNs
            Pkp1Plus(:,:,kp1,stationTracker) = eye(n);

        end % if

        % Iterate station tracker
        stationTracker = stationTracker + 1;

    end % for

    % Store perturbation calculated
    xPerturbsLKF(:,kp1) = delxkp1Plus; 

    % Store current PLKF and SLKF
    if ~isempty(visibleStationIds)
        % Calculate the combined Pkp1Minus
        PLKFMinus(:,:,kp1) = Ftilde(:,:, k) * PLKF(:,:,k) * Ftilde(:,:, k)' + Omegatilde(:,:, k) * Q * Omegatilde(:,:, k)';

        % Loop through visible stations, combine their K and Hs
        Kcombined = [];
        Hcombined = [];
        for kk = 1:length(visibleStationIds)
            Kcombined = [Kcombined, Kkp1Mat(:,:,kp1,visibleStationIds(kk))];
            Hcombined = [Hcombined; Htilde(visibleStationjjs(kk):visibleStationjjs(kk)+2, :, kp1)];
        end

        % Calculate the combined Pkp1Plus
        PLKF(:,:,kp1) = (eye(n) -  Kcombined * Hcombined) * PLKFMinus(:,:,kp1);

        % Create block diagonal big R matrix for this number of visible stations
        Rbig = [];
        for qq = 1:length(visibleStationjjs)
            Rbig = blkdiag(Rbig, R);
        end

        % Calculate SLKF
        SLKF{kp1} = ( Hcombined * PLKFMinus(:,:,kp1) * Hcombined' ) + Rbig;

        % Extract y data into one stacked vec at this timestep
        yNoisyVec = reformatYNom(yNoisy(:,kp1));

        % Extract yNom data into one stacked vec at this timestep
        yNomVec = reformatYNom(yNom(:,kp1));

        % Ensure that size of ydata here matches size of y combined
        assert(length(yNoisyVec) == size(Hcombined,1), 'ERROR: Size of ydata does not match size of yNL')

        % Calculate yError (or usual ey,k)
        yErrorLKF{kp1} = yNoisyVec - (Hcombined * delxkp1Minus + yNomVec);

        % Calculate new perturbation ey,k
        eyLKF{kp1} = (yNoisyVec - yNomVec) - (Hcombined * delxkp1Minus);

    else % No visible stations

        % Set Sk = R
        SLKF{kp1} = R;

        % Set eyk to empty
        yErrorLKF{kp1} = [];
        eyLKF{kp1} = [];

    end % if

end; clear ii jj kk qq; % for

% Calculate total states
xLKF = xNom + xPerturbsLKF;

% Calculate estimated state error using the ground truth data
xErrorLKF = xTMT - xLKF; % AKA usual ex
exLKF = (xTMT - xNom) - (xLKF - xNom); % Perturbation ex

end % function