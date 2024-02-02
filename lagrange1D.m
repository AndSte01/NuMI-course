function [f, df] = lagrange1D(points, x)
arguments
    points (:,2) double
    x (1,1) double
end
%LAGRANGE1D Evaluates a Lagrange polynomial and its derivative
% It does so by interpolating the given points and evaluating at position x
%
% Inputs:
%   points: the points to interpolate [x y; ...]
%   x:      the position x at wich the polynomial should be 
%
% Output:
%   f:  the result of the interpolated polynomial at given x
%   df: the derivative of the interpolated polynomial at given x
%
% Exercise 2.2
%
% requires lagrange_base, lagrange_base_derivative (both in this file)
%
% © 2024, Andreas Steger

% Order of Lagrange polynomial (determined by number of points)
n = size(points, 1) - 1;

% set initial values
f = 0;
df = 0;

% calculate the sum over the different base polynomials
for i = 1:(n+1)
    f = f + lagrange_base(x, n, i, i, points(:,1)) * points(i,2); % we only need to skip i, so m = i
    df = df + lagrange_base_derivative(x, n, i, points(:,1)) * points(i,2);
end
end

function val = lagrange_base_derivative(x, n, i, x_nodes)
arguments
    x (1,1) double
    n (1,1) int64
    i (1,1) int64
    x_nodes (:,1) double
end
%LAGRANGE_BASE_DERIVATIVE 
%
% Inputs:
%   x:       the position at wich the basis should be evaluated
%   n:       order of Lagrange polynomial
%   i:       first value to skip (relevant for Lagrange polynomial)
%   x_nodes: list of x values
%
% Outputs:
%   val: value of the basis polynomial
%
% © 2024, Andreas Steger

% set initial value
val = 0;

lagrange_basis_fraction_deriv = @(i,m) lagrange_base(x, n, i, m, x_nodes) / (x_nodes(i) - x_nodes(m));

for m = 1:(i-1)
    val = val + lagrange_basis_fraction_deriv(i, m);
end

for m = (i+1):(n+1)
    val = val + lagrange_basis_fraction_deriv(i, m);
end
end

function val = lagrange_base(x, n, i, m, x_nodes)
arguments
    x (1,1) double
    n (1,1) int64
    i (1,1) int64
    m (1,1) int64
    x_nodes (:,1) double
end
%LAGRANGE_BASE Evaluates the Lagrangian basis (part of lagrange polynomial)
%
% Inputs:
%   x:       the position at wich the basis should be evaluated
%   n:       order of lagrange polynomial
%   i:       first value to skip (relevant for lagrange polynomial)
%   m:       second value to skip (relevant for derivative)
%   x_nodes: list of x values
%
% Outputs:
%   val: value of the basis polynomial
%
% © 2024, Andreas Steger

% set initial value
val = 1;

for k = 1:(n+1)
    if k ~= i && k~= m
        val = val * (x - x_nodes(k)) / (x_nodes(i) - x_nodes(k));
    end
end
end