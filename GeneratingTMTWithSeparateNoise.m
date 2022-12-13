for ii = 1:nRuns
    % Generate process noise for this timestep
    wkMat = mvnrnd(zeros(length(tvec), length(x0_pert)), diag([0, Qtrue(1,1), 0, Qtrue(2,2)]));
    % Convert to column vector
    wkMat = wkMat';
    for jj = 1:nTimesteps
        xTMT(:,jj,ii) = xTMT(:,jj,ii) + wkMat(:,jj);
    end % for jj
end; clear ii; jj; % for ii