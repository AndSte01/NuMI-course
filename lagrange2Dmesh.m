function [f, df] = lagrange2Dmesh(y, xi, eta)
arguments
    y (4,1) double
    xi (:,:) double
    eta (:,:) double
end
%LAGRANGE2DMESH evaluates the lagrange2D function over a mesh of xi and eta
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
%        must be between [-1; 1]
%   eta: the eta at which to evaluate the function
%        must be between [-1; 1]
%
% Output:
%   f:  interpolated values at the given xi and eta
%   df: derivative of the interpolated function at given xi and eta
%
% © 2024, Andreas Steger

% check if sizes are correct
if(~all(size(xi) == size(eta)))
    error("xi and eta must be of same size");
end

% create empty solution arrays
dimensions = size(xi);
f = zeros(dimensions);
df = zeros([dimensions(1) dimensions(2)*2]);

% evaluate on the different positions
for k = 1:dimensions(1)
    for l = 1:dimensions(2)
        [f_temp, df_temp] = lagrange2D(y, xi(k,l), eta(k,l));
        f(k,l) = f_temp;
        df(k,2*l-1:2*l) = df_temp;
    end
end
end