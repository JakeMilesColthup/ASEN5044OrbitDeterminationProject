function yDataVec = reformatYDataCell36(ydata, nStations)
% Rearranges the current timestep's ydata cell into a stacked vector that does not include station IDs 
% Input: ydata's value of a cell array at one timestep (what is inside the cell, not the cell itself)
% Output: yDataVec which is a 36x1 vector

% Determine how many visible stations there are
nVisibleStations = size(ydata,2);

% Initialize output
yDataVec = NaN * ones(3*nStations,1);

% Loop through all stations
for ii = 1:nStations

    % Initialize count if this station has been counted
    count = 0;

    % Loop through the actually visible stations
    for jj = 1:nVisibleStations

        % Determine if this station is present in the data
        if (ydata(4,jj) == ii) && (count == 0)
            % Store this station's data
            yDataVec((ii-1)*3+1:(ii-1)*3+3) = ydata(1:3,jj);

            % Iterate count
            count = count + 1;
        end % if

    end % for jj

end; clear ii jj; % for


end % function