function [LHS,RHS] = OST(theta,timestep,M,B,C,sol)
arguments
    theta (1,1) double
    timestep (1,1) double
    M (:,:) double
    B (:,:) double
    C (:,2) double
    sol (:,1) double
end
%OST One-Step-Theta method for solving linear inhomogeneous differential
%equations.
% The equation should have the following layout:
%
%     ∂ϕ
%   M*-- = B(t)*ϕ + C(t)
%     ∂t
%
% Input:
%   theta:    Define the behavior of the equation (scalar)
%                0:   Forward-Euler-method (explicit)
%                1:   Backward-Euler-method (implicit)
%                0.5: Trapezoidal rule
%   timestep: Integration interval (scalar)
%   M:        Factor in front of dϕ/dt
%                M ∈ m×m
%   B:        Factor in front of ϕ
%                B ∈ m×m    Layout: [B(t_n+1), B(t_n)]
%   C:        Inhomogeneous part
%                C ∈ m×1    Layout: [C(t_n+1), C(t_n)]
%   sol:      solution of the timestep(s) before
%                sol ∈ m×1
%   
% Output:
%   LHS: Left-hand-side of the resulting equation system
%   RHS: Right-hand-side of the resulting equation system
%
% Exercise 6
%
% © 2024, Andreas Steger

% get height of matrices
n = height(M);

% check the dimensions
if (n ~= width(M))
    error("M must be a square matrix (is "+n+"×"+width(M)+")")
end
if (width(B) ~= 2*height(B))
    error("B must be composed of two square matrices (both currently " + ...
        "interpreted as "+height(B)+"×"+width(B)*.5+")")
end
if (n ~= height(B))
    error("M and B must be of same height (is "+height(B)+", expected: "+n+")")
end
if (height(C) ~= n)
    error("C must be composed of row vectors with the same height as M")
end
if(height(sol) ~= n)
    error("sol must be row vector with height of M")
end

% do the calculation
LHS = M - theta*timestep*B(:,1:n);
RHS = (M + (1-theta)*timestep*B(:,n+1:2*n))*sol...
    + timestep*(theta*C(:,1) + (1-theta)*C(:,2));
end