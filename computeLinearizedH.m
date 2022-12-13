function Htilde = computeLinearizedH(xstar, stationVec)
% Takes in nominal state vector and station vector information and computes 
% the corresponding nonlinear DT Jacobian H. Uses the CT Jacobians derived in problem 1.
% stationVec contains 4xn matrix where each row is Xs, Xsdot, Ys, Ysdot and each column is a new station i.

%% First, compute the CT Jacobian evaluated at the nominal values given

% Determine number of stations
nStations = size(stationVec, 2);

% Initialize size of Ctilde
Ctilde = zeros(3*nStations, 4);

% Calculate C matrix for each station and store in output
for ii = 1:nStations
    if any(isnan(stationVec(:,ii)))
        Ctilde((ii-1)*3+1:(ii-1)*3+3, :) = NaN * ones(3,4);
    else
%         Ctilde((ii-1)*3+1:(ii-1)*3+3, :) = [(xstar(1)-stationVec(1,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ), 0, (xstar(3)-stationVec(3,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ), 0;...
%           -1*( (stationVec(3,ii)-xstar(3)) * ( (stationVec(4,ii)-xstar(4)) * xstar(1) - stationVec(1,ii)*stationVec(4,ii) + (stationVec(2,ii)-xstar(2)) * stationVec(3,ii) - xstar(3)*stationVec(2,ii) + xstar(4)*stationVec(1,ii) + xstar(2)*xstar(3) ) ) / ( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 )^(3/2), ...
%           (xstar(1)-stationVec(1,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ),...
%           -1*( (stationVec(1,ii)-xstar(1)) * ( (stationVec(2,ii)-xstar(2)) * xstar(3) + (stationVec(1,ii)-xstar(1))*stationVec(4,ii) + (xstar(2)-stationVec(2,ii))*stationVec(3,ii) - xstar(4)*stationVec(1,ii) + xstar(1)*xstar(4) ) ) / ( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 )^(3/2),...
%           (xstar(3)-stationVec(3,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ); ...
%           (stationVec(3,ii)-xstar(3)) / ( (xstar(1)-stationVec(1,ii))^2 + (stationVec(3,ii)-xstar(3))^2 ), 0, (xstar(1)-stationVec(1,ii)) / ( (xstar(3)-stationVec(3,ii))^2 + (xstar(1)-stationVec(1,ii))^2 ), 0];

          Ctilde((ii-1)*3+1:(ii-1)*3+3, :) = [(xstar(1)-stationVec(1,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ), 0, (xstar(3)-stationVec(3,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ), 0;...
          
          -1*( (stationVec(3,ii)-xstar(3)) * ( (stationVec(4,ii)-xstar(4))*(xstar(1)-stationVec(1,ii)) + (stationVec(2,ii)-xstar(2))*(stationVec(3,ii)-xstar(3)) ) ) / ( (xstar(1)-stationVec(1,ii))^2 + (stationVec(3,ii)-xstar(3))^2 )^(3/2), ...
          (xstar(1)-stationVec(1,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ),...
          -1*( (stationVec(1,ii)-xstar(1)) * ( (stationVec(2,ii)-xstar(2))*(xstar(3)-stationVec(3,ii)) + (stationVec(1,ii)-xstar(1))*(stationVec(4,ii)-xstar(4)) ) ) / ( (stationVec(1,ii)-xstar(1))^2 + (xstar(3)-stationVec(3,ii))^2 )^(3/2), ...
          (xstar(3)-stationVec(3,ii)) / sqrt( (xstar(1)-stationVec(1,ii))^2 + (xstar(3)-stationVec(3,ii))^2 ); ...

          (stationVec(3,ii)-xstar(3)) / ( (xstar(1)-stationVec(1,ii))^2 + (stationVec(3,ii)-xstar(3))^2 ), 0, (xstar(1)-stationVec(1,ii)) / ( (xstar(3)-stationVec(3,ii))^2 + (xstar(1)-stationVec(1,ii))^2 ), 0];
    end
end; clear ii;

%% Convert CT Jacobian to DT

% H:
Htilde = Ctilde;


end % function