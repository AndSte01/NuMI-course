function [sysmat,rhs] = assignDBC(sysmat,rhs,dbc)
arguments
    sysmat (:,:) double
    rhs (:,1) double
    dbc (:,2) double
end
%ASSIGNDBC Function to apply dirchlet boundary condition to fem
%
% Input:
%   sysmat: The system matix
%   rhs:    The right hand side of the problem
%   dbc:    The dirichlet boundary condition
%               [node, dbc_value; ...; node, dbc_value]
%
% Output:
%   sysmat: The updated system matrix
%   rhs:    The updated right hand side
%
% Exercise 7
%
% © 2024, Andreas Steger

% iterate over all boundary conditions
for i=1:height(dbc)
    % assemble the new row
    row = zeros([1 width(sysmat)]);
    row(dbc(i,1)) = 1;

    % write row to system matrix
    sysmat(dbc(i,1),:) = row;

    % now assign the value of the boundary condition to the rhs
    rhs(dbc(i,1)) = dbc(i,2);
end

end