function [LHS,RHS] = AB2(timestep,M,B,C,sol)
arguments
    timestep (1,1) double
    M (:,:) double
    B (:,:) double
    C (:,2) double
    sol (:,2) double
end
%AB2 Adams-Bashforth-method of 2nd order for solving linear inhomogeneous
%differential equations.
% The equation should have the following layout:
%
%     ∂ϕ
%   M*-- = B(t)*ϕ + C(t)
%     ∂t
%
% Input:
%   timestep: Integration interval (scalar)
%   M:        Factor in front of dϕ/dt
%                M ∈ m×m
%   B:        Factor in front of ϕ
%                B ∈ m×m    Layout: [B(t_n), B(t_n-1)]
%   C:        Inhomogeneous part
%                C ∈ m×1    Layout: [C(t_n), C(t_n-1)]
%   sol:      solution of the timestep(s) before
%                sol ∈ m×1  Layout: [ϕ(t_n), ϕ(t_n-1)]
%
% Output:
%   LHS: Left-hand-side of the resulting equation system
%   RHS: Right-hand-side of the resulting equation system
%
% Exercise 6
%
% © 2024, Andreas Steger

% get width of matrices
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
    error("M and B must be of same height ("+n+" and "+height(B)+")")
end
if(height(C) ~= n)
    error("C must be composed of row vectors with height of M")
end
if(height(sol) ~= n)
    error("sol must be composed of row vectors with height of M")
end

% do the calculation
LHS = M;
RHS = M*sol(:,1) + timestep/2*(...
    3*( B(:,1:n)*sol(:,1) + C(:,1) )-...
    ( B(:,n+1:2*n)*sol(:,2) + C(:,2)));
end