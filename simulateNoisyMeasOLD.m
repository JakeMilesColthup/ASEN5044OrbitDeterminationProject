function yNoisy = simulateNoisyMeasOLD(ydata, R, nStations)
% Simulates noisy measurements by adding noise to the provided ydata according to given R matrix from Canvas.
% Simulates for all timesteps, but only for visible stations (where data already exists in those elements)

% Calculate Sv
Sv = chol(R, 'lower');

% % Set seed to replicate results temporarily
% rng(100);

% Initialize output
yNoisy = NaN * ones(3*nStations, length(ydata));

% Simulate noisy measurements for each timestep

% Loop through each timestep
for ii = 1:length(ydata)
    % Determine how many stations have data at this timestep
    nVisibleStations = length(ydata{ii}(~isnan(ydata{ii}))) / 4;

    % Loop through each visible station
    for jj = 1:nVisibleStations
        stationId = ydata{ii}(4,jj);
        vk = Sv * randn(3,1);
        yNoisy((stationId*3) - 2 : (stationId*3), ii) = ydata{ii}(1:3,jj) + vk;
    end 
    
end; clear ii jj; % for

end % function