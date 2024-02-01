function [f, df] = lagrange2D(y, xi, eta)
arguments
    y (4,1) double
    xi (1,1) double
    eta (1,1) double
end
%LAGRANGE2D Evaluates a 2D Lagrange polynomial and its derivative
% The function works on the given definition of i and points
% 
%  [-1,1]     [1,1]
%     4 ------- 3
%     |         |
%     |         |
%     |         |
%     1 --------2
%  [-1,-1]    [1,-1]
%
% Inputs:
%   y:   the function values of f at the corners in the order stated above
%   xi:  the xi at which to evaluate the function
%        must be between [-1; 1], must be a scalar
%   eta: the eta at which to evaluate the function
%        must be between [-1; 1], must be a scalar
%
% Outputs:
%   f:  interpolated function values
%   df: derivatives at the given point [dxi, deta]
%
% requires linquadref linquadderivref
%
% © 2024, Andreas Steger

%% bring input to form
% make sure y is always a column vector
% if isrow(y)
%     y = y';
% end

%size(xi)
%size(eta)

%% evaluate f
f = linquadref(xi, eta)' * y;

%% calculate df
% empty variable for df
df = zeros([1 2]);

% calculate the base polynomials
base_poly_deriv = linquadderivref(xi, eta);

% evaluate df
df(:,1) = base_poly_deriv(:,1)' * y;
df(:,2) = base_poly_deriv(:,2)' * y;
end