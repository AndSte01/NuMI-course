function [equal] = feq(a, b, varargin)
%FEQ Compares two floating point numbers and determines if they are equal
% a and b can be vectors (maybe even matrices) of same size, evaluation is
% then element wise
%
% Inputs:
%   a: first number
%   b: second number
%   d: max difference (optional)
%
% © 2024, Andreas Steger

% check if the user provided to many input arguments
if nargin > 3
    error("Too many input parameters (max. 3)");
end

% decide wether to use the embedded delta or the one provided by the user
% difference_max: maximum difference between a, b up to which both are considered equal
if nargin > 2
    % get the delta from user input
    difference_max = varargin{1};
else
    % constants for configuration
    difference_max = 1e-12;
end

% logic
equal = all(abs(a-b) < difference_max, 'all');
end
