function [stationVec, stationIDs] = checkVisibleStationsNoisyMeas(nStations, RE, omegaE, t, yNoisy)
% Checks which stations can see the satellite at this timestep, computes vector of station parameters for each of those stations
% Input: x is current state vector at this timestep, yNoisy at this timestep (36x1 with NaNs)
% Output: 4xn where n is number of stations that can see the satellite

% Initialize outputs
stationVec = NaN * ones(4, nStations); % column has Xs, Xsdot, Ys, Ysdot
stationIDs = []; % Reports station IDs of visible stations

% Loop through and calculate stationVec
for ii = 1:nStations
    % Calculate this station's parameters
    thetai0 = (ii-1) * (pi/6);
    Xi = RE * cos(omegaE*t + thetai0);
    Yi = RE * sin(omegaE*t + thetai0);
    Xdoti = omegaE * -RE * sin(omegaE*t + thetai0);
    Ydoti = omegaE * RE * cos(omegaE*t + thetai0);

    % Check if this station can see the satellite
    if ~isnan(yNoisy(3*ii-2))
        % It can see, so add it to the nominal stationVec
        stationVec(:,ii) = [Xi; Xdoti; Yi; Ydoti];

        % Also store its id
        stationIDs = [stationIDs, ii];
    end % if

end; clear ii; % for

end % function