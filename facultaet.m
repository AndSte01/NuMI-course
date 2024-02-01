function nfac = facultaet(n)
arguments
    n int64
end
%facultaet Calculates the faculty of a given number
%
% Inputs:
%   n: The number to calculate the faculty from
%
% Exercise 1.2
%
% © 2024, Andreas Steger

if n ~= 0
    nfac = n * facultaet(n - 1);
else
    nfac = 1;
end
end