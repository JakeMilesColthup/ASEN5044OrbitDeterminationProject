function [stationVec, stationIDs] = checkVisibleStations(nStations, RE, omegaE, x, t)
% Checks which stations can see the satellite at this timestep, computes vector of station parameters for each of those stations
% Input: x is current state vector at this timestep
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

    % Calculate thetait
    thetait = atan2(Yi,Xi);

    % Calculate phii
    phii = atan2( (x(3) - Yi),(x(1) - Xi) );

    % Ensure that the proper limits are constructed to match the sign of phi
    if phii < 0 % Phi is negative
        if (pi/2 + thetait > pi)
            % Phi is negative, but the upper limit is too positive, so limits must be wrapped around to be negative

            % Define upper and lower limits
            lowerLim = -pi/2 + thetait - 2*pi;
            upperLim = pi/2 + thetait - 2*pi;
        else
            % Otherwise, upper limit is not too positive, so limits can be defined as normal

            % Define upper and lower limits
            lowerLim = -pi/2 + thetait;
            upperLim = pi/2 + thetait;
        end % if

    else % Phi is positive
        if (-pi/2 + thetait < - pi)
            % Phi is positive, but the lower limit is too negative, so limits must be wrapped around to be postive

            % Define upper and lower limits
            lowerLim = -pi/2 + thetait + 2*pi;
            upperLim = pi/2 + thetait + 2*pi;
        else
            % Otherwise, lower limit is not too negative, so limits can be defined as normal

            % Define upper and lower limits
            lowerLim = -pi/2 + thetait;
            upperLim = pi/2 + thetait;
        end % if
    end % if

    % Check if this station can see the satellite
    if (phii >= lowerLim) && (phii <= upperLim)
        % It can see, so add it to the nominal stationVec
        stationVec(:,ii) = [Xi; Xdoti; Yi; Ydoti];

        % Also store its id
        stationIDs = [stationIDs, ii];
    end % if

end; clear ii; % for

end % function