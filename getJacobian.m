function [J, detJ, invJ] = getJacobian(nodes, xi, eta)
arguments
    nodes (:,2) double
    xi (1,:) double
    eta (1,:) double
end
%GETJACOBIAN Calculates the Jacobian Matrix of the isoparametrical
%transformation of xi, eta to x (see getxPos) for details
% This function uses an numerical solution
%
% Inputs:
%   nodes: Nodes list determining the desired transformation space
%              layout: [nodes(1,1), nodes(1,2); ...; nodes(4,1), nodes(4,2)]
%   xi:    First coordinate in reference space
%   eta:   Second coordinate in reference space
%
% Output:
%   J:    Jacobian
%   detJ: Determinant of Jacobian
%   invJ: Inverse of Jacobian
%
% Exercise 5
%
% requires linquadderivref
%
% © 2024, Andreas Steger

% some error checking
if (length(xi) ~= length(eta))
    error("lengths of xi and eta must match");
end
if (height(nodes) ~= 4*length(xi)) % number of elements in linquadderivref
    error("not enough nodes for given xi and eta (expected "+4*length(xi)+ ...
        " got "+height(nodes)+")");
end

% number of results to calculate
n = length(xi);

% do calculations
J = zeros([2*n 2]);
detJ = zeros([n 1]);
invJ = zeros([2*n 2]);

for i = 1:n
    % calculate the Jacobian for the distinct pair of xi, eta
    J_tmp = nodes' * linquadderivref(xi, eta);
    % store result in result matrix
    J((2*n-1):(2*n),:) = J_tmp;

    % using built in function (using exact analytical solution would be to
    % expensive (and would get to much on my nerves!!!)
    detJ(i) = det(J_tmp);
    invJ((2*n-1):(2*n),:) = inv(J_tmp);
end
end