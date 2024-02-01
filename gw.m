function gaussw = gw(n)
arguments
    n int64
end
%GW Weights for 1D-Gauß-integration in form of row-vector
%
% Input:
%   n: Number of points
%
% Output:
%   gaussw: weights for integration
%
% Exercise 5
%
% requires gxw
%
% © 2024, Andreas Steger

[~, gaussw] = gxw(n);
end