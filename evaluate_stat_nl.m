function [elemat,elevec] = evaluate_stat_nl(elenodes,gpx,gpw,elesol,q,dq,optional)
arguments
    elenodes (4,2) double
    gpx (:,2) double
    gpw (:,1) double
    elesol (:,1) double
    q
    dq
    optional.lambda double = 48
end
%EVALUATE_STAT Calculates the values of the element matrices of the given
%element.
% Solving a certain stationary thermal conduction problem
%
%   -∇⋅(λ∇T) = q̇(T)
%
% Input:
%   elenode:       List of points creating the element
%                      [x, y; ...; x, y]
%   gpx:           Positions ξ_i for the gauß-integration
%   gpw:           Weights w_i for the gauß-integration
%
% Optional: (optional.x)
%   lambda:       Thermal conductivity of the material [W/mK]
%
% Output:
%   elemat: Elementmatrix
%   elevec: Elementvector
%
% Exercise 10.2
%
% requires getJacobian, linquadref, linquadderivref
%
% © 2024, Andreas Steger

% get the number of iteration steps
n_k = length(gpw);
n_N = 4; % length of output from linquadderivref

%% calculate some values that we then use to assemble the elevec and elemat
% placeholders for the temporary values
temp_v1 = zeros([n_N 1]);
temp_M1 = zeros([n_N, n_N]);
temp_M2 = zeros([n_N, n_N]);

% iterate over the gauss points
for k=1:n_k
    % [gpx(k,1), gpx(k,2)]
    % evaluate the Lagrange polynoms at the given position
    dN = linquadderivref(gpx(k,1), gpx(k,2));
    N = linquadref(gpx(k,1), gpx(k,2));
    [~, detJ, invJ] = getJacobian(elenodes, gpx(k,1), gpx(k,2));
    % iterate over rows
    for i=1:n_N
        % Calculate the value for usage in q or dq
        temp_sum_j_v1 = N'*elesol;
        % iterate over columns
        for j=1:n_N
            % N(i,:)*invJ: *invJ is important since we evaluate the
            % Lagrange function on the reference element but we need to
            % calculate it on the real element (nachdifferenzieren)
            temp_M1(i,j) = temp_M1(i,j) + (dN(i,:)*invJ*(dN(j,:)*invJ)')*detJ*gpw(k);
            % second matrix
            temp_M2(i,j) = temp_M2(i,j) + N(i,:)*N(j,:)*dq(temp_sum_j_v1)*detJ*gpw(k);
        end
        temp_v1(i) = temp_v1(i) + N(i)*q(temp_sum_j_v1)'*detJ*gpw(k);
    end
end

% apply lambda to first element
temp_M1 = optional.lambda*temp_M1;


%% generate elevec and elemat
elevec = temp_M1*elesol - temp_v1;
elemat = temp_M1 - temp_M2;

end