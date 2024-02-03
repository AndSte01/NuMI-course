function [elemat, elevec] = evaluate_instat(elenodes, gpx, gpw, elesol, eleosol, timInt_m, timestep, theta, firststep, optional)
arguments
    elenodes (4,2) double
    gpx (:,2) double
    gpw (:,1) double
    elesol (:,1) double
    eleosol (:,1) double
    timInt_m (1,1) int64 {mustBeMember(timInt_m,[1,2,3,4])}
    timestep (1,1) double
    theta (1,1) double
    firststep
    optional.lambda (1,1) double = 48
    optional.rho (1,1) double = 7800
    optional.c (1,1) double = 452
    optional.q_dot (1,1) double = 0
end
%EVALUATE_INSTAT Calculates the values of the element matrices of the given
%element.
% Solving a certain instationary thermal conduction problem
%
%   ρc(∂T/∂t) - ∇⋅(λ∇T) = q̇
%
% Input:
%   elenode:   List of points creating the element
%                  [x, y; ...; x, y]
%   gpx:       Positions ξ_i for the gauß-integration
%   gpw:       Weights w_i for the gauß-integration
%   elesol:    Solution at time (n)
%   eleosol:   Solution at time (n-1)
%   timInt_m:  Selection of the time integration method
%                  1 = OST, 2 = AB2, 3 = AM3, 4 = BDF2
%   timestep:  Timestep for time integration
%   theta:     Theta for one-step-theta integration
%                  (only relevant if timInt_m = 1)
%   firststep: Not implemented yet
%
% Optional: (optional.x)
%   lambda: Thermal conductivity of the material [W/mK]
%   rho:    Density of the material [kg/m³]
%   c:      Specific thermal capacity of the material [J/kgJ]
%   q_dot:  Evenly distributed stationary thermal source
%
% Output:
%   elemat: Elementmatrix
%   elevec: Elementvector
%
% Exercise 8
%
% requires getJacobian (linquadderivref), linquadderivref
%
% © 2024, Andreas Steger

% get the number of gauss points
n_k = length(gpw);
n_N = 4; % length of output from linquadderivref

% generate inhomogenus part
C = optional.q_dot * ones([n_N 1]);

% empty result matrix
M = zeros([n_N n_N]);
B = zeros([n_N n_N]); % in stationary state we don't set the element vector to zero

% iterate over the gauss points
for k = 1:n_k
    % [gpx(k,1), gpx(k,2)]
    % evaluate the Lagrange polynoms (and their derivatives) at the given
    % location
    N = linquadref(gpx(k,1), gpx(k,2));
    dN = linquadderivref(gpx(k,1), gpx(k,2));
    [~,detJ,invJ] = getJacobian(elenodes, gpx(k,1), gpx(k,2));

    % Build the matrix
    for i = 1:n_N
        for j = 1:n_N
            % Calculate the M (for time integration)
            M(i,j) = M(i,j) + (N(i) * N(j) * detJ) * gpw(k);

            % N(i,:)*invJ: *invJ is important since we evaluate the
            % Lagrange function on the reference element but we need to
            % calculate it on the real element (nachdifferenzieren)

            % Calculate the B (for time integration)
            B(i,j) = B(i,j) - ((dN(i,:)*invJ * (dN(j,:)*invJ)') * detJ) * gpw(k);
        end
    end
end

% multiply prefactors for M
M = optional.c * optional.rho * M;

% multiply prefactors for B
B = optional.lambda * B;

% Do the time integration
switch timInt_m
    case 1 % OST
        [elemat,elevec] = OST(theta, timestep, M, [B B], [C C], elesol);

    case 2 % AB2
        [elemat,elevec] = AB2(timestep, M, [B B], [C C], [elesol eleosol]);

    case 3 % AM3
        [elemat,elevec] = AM3(timestep, M, [B B B], [C C C], [elesol eleosol]);

    case 4 % BDF2
        [elemat,elevec] = BDF2(timestep, M, B, C, [elesol eleosol]);

    otherwise
        error("%d is not a viable option of time integration methods\n" + ...
            "    Possible options are: 1 = OST, 2 = AB2, 3 = AM3, 4 = BDF2", timInt_m);
end
end