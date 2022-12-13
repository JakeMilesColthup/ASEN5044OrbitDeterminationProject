function yNomVec = reformatYNom(yNom)
% Rearranges the current timestep's yNom NaN vector into a stacked vector without the NaNs 
% Input: yNom 36x1 vector at one timsetep
% Output: yNomVec which is a stack of 3x1 y vectors

% Initialize output
yNomVec = [];

% Store values that aren't NaN
for ii = 1:length(yNom)
    % Determine if this element is NaN or not
    if ~isnan(yNom(ii))

        % Store this data
        yNomVec = [yNomVec; yNom(ii)];

    end %if
end; clear ii; % for

end % function