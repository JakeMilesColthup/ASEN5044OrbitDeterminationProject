function yDataOneStation = reformatYDataCellOneStation(ydata,stationID)
% Rearranges the current timestep's ydata cell into a 3x1 vector that doesn't include stationID
% Input: ydata's value of a cell array at one timestep (what is inside the cell, not the cell itself)
% Output: yDataOneStation which is a of 3x1 vector


% Loop through and find the column corresponding to this stationID
for ii = 1:size(ydata,2)
    if ydata(4,ii) == stationID
        % Store this station's ydata
        yDataOneStation = ydata(1:3,ii);
    end % if
end; clear ii; % for

assert(~isempty(yDataOneStation), 'Station ID not found in data')

end % function