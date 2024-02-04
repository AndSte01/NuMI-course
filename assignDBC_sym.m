function [sysmat,rhs] = assignDBC_sym(sysmat,rhs,dbc)
arguments
    sysmat (:,:) double
    rhs (:,1) double
    dbc (:,2) double
end
%ASSIGNDBC Function to apply dirchlet boundary condition to fem while
%keeping symmetric properties of system matrix.
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
% Exercise 9
%
% © 2024, Andreas Steger

% iterate over all boundary conditions
for i=1:height(dbc)
    % assign boundary condition value to all rows
    rhs = rhs - dbc(i,2)*sysmat(:,dbc(i,1));

    % set boundary condition value on selected row
    rhs(dbc(i,1)) = dbc(i,2);

    % assemble the new row
    row = zeros([1 width(sysmat)]);
    row(dbc(i,1)) = 1;
    
    % set column to zeros and write row to system matrix
    sysmat(:,dbc(i,1)) = 0;
    sysmat(dbc(i,1),:) = row;
end
end