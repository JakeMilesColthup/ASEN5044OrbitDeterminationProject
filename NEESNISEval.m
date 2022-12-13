function NEESNISEval(alpha, r1x, r2x, r1y, r2y, epsilonxBar, epsilonyBar)
% Test to see if 1-alpha% within NIS/NEES

% Define goal percent
goal = 1-alpha; 

% Loop through and count how many NEES stats lie in the ranges
NEEScount = 0;
for ii = 1:length(r1x)
    if (epsilonxBar(ii) > r1x(ii)) && (epsilonxBar(ii) < r2x(ii))
        NEEScount = NEEScount + 1;
    end
end; clear ii; % for

% Loop through and count how many NIS stats lie in the ranges
NIScount = 0;
for ii = 1:length(r1y)
    if (epsilonyBar(ii) > r1y(ii)) && (epsilonyBar(ii) < r2y(ii))
        NIScount = NIScount + 1;
    end
end; clear ii; % for

% Calculate if results meet expectations
NEESsuccess = NEEScount/length(r1x) > goal;
NISsuccess = NIScount/length(r1y) > goal;

% Print results
if NEESsuccess
    fprintf('NEES Passed! \n')
else
    fprintf('NEES Failed. \n')
end

if NISsuccess
    fprintf('NIS Passed! \n')
else
    fprintf('NIS Failed. \n')
end

end % function