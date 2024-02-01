function gaussw = gw2dref(n)
arguments
    n int64
end
%GXW2DREF Positions weights for 2D-Gauß-integration
%
% Input:
%   n: Number of points
%
% Output:
%   gaussw: weights for 2D-integration
%
% Exercise 5
%
% requires gxw2dref  (gxw)
%
% © 2024, Andreas Steger

% only return desired value
[~, gaussw] = gxw2dref(n);
end