function yNoisy = simulateNoisyMeas(yNL, R, nStations)
% Simulates noisy measurements by adding noise to the nonlinear measurement model according to given R matrix from Canvas.
% Simulates for all timesteps, but only for visible stations (where data already exists in those elements)

% Calculate Sv
Sv = chol(R, 'lower');

% % Set seed to replicate results temporarily
% rng(100);

% Initialize output
yNoisy = NaN * ones(3*nStations, length(yNL));

% Simulate noisy measurements for each timestep

% Loop through each timestep
for ii = 1:length(yNL)

    % Loop through each station
    for jj = 1:3:3*nStations

        % Check if this station is visible
        if all(~isnan(yNL(jj:jj+2, ii)))
            % It is visible, so generate noisy measurements
            vk = Sv * randn(3,1);
            yNoisy(jj:jj+2, ii) = yNL(jj:jj+2, ii) + vk;
        end % if

    end % for
    
end; clear ii jj; % for

end % function