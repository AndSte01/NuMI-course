function deriv = linquadderivref(xi,eta)
arguments
    xi (1,:) double
    eta (1,:) double
end
%LINQUADDERIVREF Evaluates Lagrangian base polynomials for given xi and eta
% The function works on the given definition of i and points
% 
%  [-1,1]    [1,1]
%     4 ------ 3
%     |        |
%     |        |
%     |        |
%     1 -------2
%  [-1,-1]   [1,-1]
%
% Input:
%   xi:  the xi at which to evaluate (must be between [-1; 1])
%   eta: the eta at which to evaluate (must be between [-1; 1]
%
% Output:
%   val: evaluated base polynomials in the order [N_1_dxi, N_1_deta; ...; N_4_dxi, N_4_deta]
%        in counter clockwise enumeration starting at the bottom left
%        corner. If xi and eta are vectors the first half of the columns
%        will contain all the dxi's and the second half all the deta's.
%
% Exercise 3
%
% © 2024, Andreas Steger

% derivatives of the base polynomials
N_i_dxi = @(eta, xi_0, eta_0) 0.25 * (eta + eta_0) * sign(xi_0 * eta_0);
N_i_deta = @(xi, xi_0, eta_0) 0.25 * (xi + xi_0) * sign(xi_0 * eta_0);

% evaluate the derivatives of the Lagrangian base polynomials
deriv = [ N_i_dxi(eta, -1, -1) N_i_deta(xi, -1, -1);
    N_i_dxi(eta, 1, -1) N_i_deta(xi, 1, -1);
    N_i_dxi(eta, 1, 1) N_i_deta(xi, 1, 1);
    N_i_dxi(eta, -1, 1) N_i_deta(xi, -1, 1);];
end