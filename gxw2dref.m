function [gaussx, gaussw] = gxw2dref(n)
arguments
    n int64
end
%GXW2DREF Positions eta_i and weights for 2D-Gauß-integration
%
% Input:
%   n: Number of points
%
% Output:
%   gaussx: xi for 2D-integration
%   gaussw: weights for 2D-integration
%
% Exercise 5
%
% requires gxw
%
% © 2024, Andreas Steger

% get weights and etas in 1D
[x, w] = gxw(n);

% combine the wights and etas to different points
gaussx = table2array(combinations(x, x));
gaussw = prod(table2array(combinations(w, w)), 2);

end