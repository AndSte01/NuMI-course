function [gaussx, gaussw] = gxw(n)
arguments
    n int64
end
%GXW Positions xi_i and weights for 1D-Gauß-integration in
%form of row-vector
%
% Input:
%   n: Number of points
%
% Output:
%   gaussx: xi for integration
%   gaussw: weights for integration
%
% Exercise 5
%
% © 2024, Andreas Steger

% Look up table for the different numbers of points
switch n
    case 1
        gaussx = 0;
        gaussw = 2;

    case 2
        gaussx = [-1/sqrt(3) 1/sqrt(3)];
        gaussw = [1 1];

    case 3
        gaussx = [-sqrt(3/5) 0 sqrt(3/5)];
        gaussw = [5/9 8/9 5/9];

    case 4
        gaussx = [-sqrt(3/7 + 2/7*sqrt(6/5)), -sqrt(3/7 - 2/7*sqrt(6/5)),...
            sqrt(3/7 - 2/7*sqrt(6/5)), sqrt(3/7 + 2/7*sqrt(6/5))];
        gaussw = [(18 - sqrt(30))/36, (18 + sqrt(30))/36,...
            (18 + sqrt(30))/36, (18 - sqrt(30))/36];

    case 5
        gaussx = [-1/3*sqrt(5+2*sqrt(10/7)), -1/3*sqrt(5-2*sqrt(10/7)),...
            0, 1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];
        gaussw = [(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225,...
            (322+13*sqrt(70))/900, (322-13*sqrt(70))/900];

    otherwise
        error("Only integers between 1 and 5 are supported");
end
end