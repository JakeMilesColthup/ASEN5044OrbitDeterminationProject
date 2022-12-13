function [Ftilde, Gtilde, Omegatilde] = computeLinearizedDyn(xstar, deltaT, mu)
% Takes in nominal state vector, control inputs, and deltaT and computes the DT Linearized Model Matrices from 
% the corresponding nonlinear DT Jacobians. Uses the CT Jacobians derived in problem 1.

%% First, compute the CT Jacobians evaluated at the nominal values given

% A:
Atilde = [0 1 0 0; ...
         (mu*(2*xstar(1)^2 - xstar(3)^2))/(xstar(1)^2 + xstar(3)^2)^(5/2), 0, (3*mu*xstar(1)*xstar(3))/(xstar(1)^2 + xstar(3)^2)^(5/2), 0;...
         0, 0, 0, 1; ...
         (3*mu*xstar(1)*xstar(3))/(xstar(1)^2 + xstar(3)^2)^(5/2), 0, (mu*(2*xstar(3)^2 - xstar(1)^2))/(xstar(1)^2 + xstar(3)^2)^(5/2), 0];

% B:
Btilde = [0 0; 1 0; 0 0; 0 1];

% Gamma:
Gammatilde = [0 0; 1 0; 0 0; 0 1];

%% Convert CT Jacobians to DT

% F:
Ftilde = eye(size(Atilde)) + deltaT * Atilde;

% G:
Gtilde = deltaT * Btilde;

% % Create Ahat
% Ahat = [Atilde Btilde; zeros(2,6)];
% 
% % Perform matrix exponential calculation
% exponentialMat = expm(Ahat * deltaT);
% 
% % Extract and display the F and G matrices (rounded to 4 decimals)
% Ftilde = exponentialMat(1:size(Atilde,1), 1:size(Atilde,2));
% Gtilde = exponentialMat(1:size(Btilde,1), size(Atilde,2)+1:end);

% Omega:
Omegatilde = deltaT * Gammatilde;

end % function