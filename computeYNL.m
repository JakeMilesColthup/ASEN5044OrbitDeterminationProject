function yNL = computeYNL(nStations, stationIDs, RE, omegaE, x, t)
% Checks which stations can see the satellite at this timestep, computes vector of station parameters for each of those stations
% Input: x is current state vector at this timestep, stationIDs is which stations are visible at this timestep
% Output: yNL at this timestep (36x1) with potential NaNs if station isn't visible

% Initialize outputs
yNL = NaN * ones(3*nStations, 1);

% Loop through and calculate stationVec
for ii = 1:nStations
    % Calculate this station's parameters
    thetai0 = (ii-1) * (pi/6);
    Xi = RE * cos(omegaE*t + thetai0);
    Yi = RE * sin(omegaE*t + thetai0);
    Xdoti = omegaE * -RE * sin(omegaE*t + thetai0);
    Ydoti = omegaE * RE * cos(omegaE*t + thetai0);

    % Calculate phii
    phii = atan2( (x(3) - Yi),(x(1) - Xi) );

    % Check if this station can see the satellite
    if any(stationIDs == ii)
        % Compute rho
        rhoi = sqrt( (x(1)-Xi)^2 + (x(3)-Yi)^2 ); 

        % Compute rho dot
        rhodoti = ( (x(1)-Xi)*(x(2)-Xdoti) + (x(3)-Yi)*(x(4)-Ydoti) ) / rhoi;

        % Store yNL
        yNL((ii*3) - 2 : (ii*3)) = [rhoi; rhodoti; phii];
    end % if

end; clear ii; % for

end % function