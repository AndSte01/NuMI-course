function [sysmat,rhs] = assignDBC_nl(sysmat,rhs,dbc)
arguments
    sysmat (:,:) double
    rhs (:,1) double
    dbc (:,2) double
end
%ASSIGNDBC Function to apply dirichlet boundary condition to stationary
%inhomogeneous problem
%
% Input:
%   sysmat: The system matrix
%   rhs:    The right hand side of the problem
%   dbc:    The dirichlet boundary condition
%               [node_1, ..., node_n]
%
% Output:
%   sysmat: The updated system matrix
%   rhs:    The updated right hand side
%
% Helper for Exercise 10
%
% © 2024, Andreas Steger

% iterate over all boundary conditions
for i=1:height(dbc)
    % assemble the new row
    row = zeros([1 width(sysmat)]);
    row(dbc(i,1)) = 1;

    % write row to system matrix
    sysmat(dbc(i,1),:) = row;

    % set row in the rhs to zero
    rhs(dbc(i,1)) = 0;
end
end