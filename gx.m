function gaussx = gx(n)
arguments
    n int64
end
%GX Positions xi_i for 1D-Gauß-integration in form of row-vector
%
% Input:
%   n: Number of points
%
% Output:
%   gaussx: xi for integration
%
% Exercise 5
%
% requires gxw
%
% © 2024, Andreas Steger

[gaussx, ~] = gxw(n);
end