function [elemat,elevec] = evaluate_stat(elenodes,gpx,gpw,optional)
arguments
    elenodes (4,2) double
    gpx (:,2) double
    gpw (:,1) double
    optional.lambda (1,1) double = 48
end
%EVALUATE_STAT Calculates the values of the element matrices of the given
%element.
% Solving a certain stationary thermal conduction problem
%
%   -∇⋅(λ∇T) = 0
%
% Input:
%   elenode:       List of points creating the element
%                      [x, y; ...; x, y]
%   gpx:           Positions ξ_i for the gauß-integration
%   gpw:           Weights w_i for the gauß-integration
%
% Optional:
%   lambda:       Thermal conductivity of the material [W/mK]
%
% Output:
%   elemat: Elementmatrix
%   elevec: Elementvector
%
%
% Exercise 7
%
% requires getJacobian (linquadderivref), linquadderivref
%
% © 2024, Andreas Steger

% get the number of iteration steps
n_k = length(gpw);
n_N = 4; % length of output from linquadderivref

% empty result matrix
elemat = zeros([n_N, n_N]);
elevec = zeros([n_N 1]); % in stationary state we don't set the element vector to zero

% iterate over the gauss points
for k=1:n_k
    % [gpx(k,1), gpx(k,2)]
    % evaluate the Lagrange polynoms at the given position
    N = linquadderivref(gpx(k,1), gpx(k,2));
    [~, detJ, invJ] = getJacobian(elenodes, gpx(k,1), gpx(k,2));
    for i=1:n_N
        for j=1:n_N
            % N(i,:)*invJ: *invJ is important since we evaluate the
            % Lagrange function on the reference element but we need to
            % calculate it on the real element (nachdifferenzieren)
            elemat(i,j) = elemat(i,j) + (N(i,:)*invJ*(N(j,:)*invJ)')*detJ*gpw(k);
        end
    end
end

% add lambda
elemat = optional.lambda*elemat;
end