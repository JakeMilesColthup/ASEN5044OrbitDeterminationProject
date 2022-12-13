function yDataVec = reformatYDataCell(ydata)
% Rearranges the current timestep's ydata cell into a stacked vector that does not include station IDs 
% Input: ydata's value of a cell array at one timestep (what is inside the cell, not the cell itself)
% Output: yDataVec which is a stack of 3x1 y vectors

% Determine how many visible stations there are
nVisibleStations = size(ydata,2);

% If no visible stations, return a NaN 3x1 vector
if all(isnan(ydata)) 
    yDataVec = NaN * ones(3,1);
else
    % If not completely empty, loop through each visible station and create yDataVec
    yDataVec = [];

    for ii = 1:nVisibleStations
        % Store this station's ydata
        yDataVec = [yDataVec; ydata(1:3,ii)];
    end; clear ii; % for

end % if


end % function