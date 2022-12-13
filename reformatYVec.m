function yOneStation = reformatYVec(yVec, visibleStationjj)
% Uses the current timestep's stacked yVec to create a new vector of just the value for that station
% Input: yVec 36x1 vector at one timsetep
% Output: yOneStation which is a 3x1 vector

yOneStation = yVec(visibleStationjj:visibleStationjj+2);

end % function