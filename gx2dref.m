function gaussx = gx2dref(n)
arguments
    n int64
end
%GX2DREF Positions eta_i for 2D-Gauß-integration
%
% Input:
%   n: Number of points
%
% Output:
%   gaussx: xi for 2D-integration
%
% Exercise 5
%
% requires gxw2dref (gxw)
%
% © 2024, Andreas Steger

% only return desired value
gaussx = gxw2dref(n);
end