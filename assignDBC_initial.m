function T = assignDBC_initial(T,dbc)
arguments
    T (:,1) double
    dbc (:,2) double
end
%ASSIGNDBC_INITIALa assign dbc to initial conditions
%
% Input:
%   T:   Initial values
%   dbc: The dirichlet boundary condition
%            [node_1, ..., node_n]
%
% Output:
%   T: Initial values fulfilling the boundary condition
%
% Helper for Exercise 10
%
% © 2024, Andreas Steger

% iterate over all boundary conditions
for i=1:height(dbc)
    % now assign the value of the boundary condition to the rhs
    T(dbc(i,1)) = dbc(i,2);
end
end