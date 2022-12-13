function [yNomVec, yNoisyVec] = reformatYNomWithNoisy(yNom, yNoisy)
% Rearranges the current timestep's yNom 36x1 vector and the current timestep's yNoisy 36x1 vector with NaNs
% Input: yNom 36x1 vector at one timsetep
% Output: yNomVec which is a stack of 3x1 y vectors

% Initialize output
yNomVec = [];
yNoisyVec = [];

% Store values that aren't NaN
for ii = 1:length(yNoisy)
    % Determine if this element is NaN or not
    if ~isnan(yNoisy(ii))

        % Store this data
        yNomVec = [yNomVec; yNom(ii)];
        yNoisyVec = [yNoisyVec; yNoisy(ii)];

    end %if
end; clear ii; % for

assert(length(yNomVec) == length(yNoisyVec), 'Size of yNomVec and yNoisyVec are different')

end % function