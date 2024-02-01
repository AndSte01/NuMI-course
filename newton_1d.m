function [x, n] = newton_1d(f, df, x0, dx, n)
arguments
    f
    df
    x0 (1,1) double
    dx (1,1) double
    n (1,1) int64
end
%NEWTON_1D Newtons method for iteratively finding roots of functions
%
% Input:
%   f:  function f of which to find the root
%   df: first derivative of f
%       (might be implemented as finite difference scheme)
%   x0: start point of iterations
%   dx: minimal difference between two iteration steps
%       (accuracy of the method)
%   n:  Maximum number of iterations to perform
%
% Output:
%   x:  the root or the value after the last iteration
%
% Exercise 10.1
%
% © 2024, Andreas Steger

% set start value
x = x0;

% now do the iterations till we reach the maximum or the error was small
% enough
for i=1:n
    x_pref = x;
    x = x - f(x)/df(x);
    
    % stop if error gets small
    if abs(x_pref - x) < dx 
        break
    end
end
end