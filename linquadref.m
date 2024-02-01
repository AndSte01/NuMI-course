function val = linquadref(xi,eta)
arguments
    xi (1,:) double
    eta (1,:) double
end
%LINQUADREF Evaluates Lagrangian base polynomials for given xi and eta
% The function works on the given points
% 
%  [-1,1]     [1,1]
%     4 ------- 3
%     |         |
%     |         |
%     |         |
%     1 --------2
%  [-1,-1]    [1,-1]
%
% Input:
%   xi:  the xi at which to evaluate (must be between [-1; 1])
%   eta: the eta at which to evaluate (must be between [-1; 1]
%
% Output:
%   val: evaluated base polynomials in the order [N_1; N_2; N_3; N_4] (as
%        column vector) in counter clockwise enumeration starting at the
%        bottom left corner
%
% Exercise 3
%
% © 2024, Andreas Steger

% if Matlab would support macros it would probably increase the
% performance to implement the N_i function as such

% defining the base polynomials
% constraints:
%   - abs(xi_0) and abs(eta_0) must be 1 (else constant 0.25 isn't correct)
%   - the last part of the product is obviously for the correct sign
N_i = @(xi, eta, xi_0, eta_0) 0.25 * (xi + xi_0) .* (eta + eta_0) * sign(xi_0 * eta_0);

% return the different base polynomials from the general description
val = [N_i(xi, eta, -1, -1); N_i(xi, eta, 1, -1); N_i(xi, eta, 1, 1); N_i(xi, eta, -1, 1)];
end