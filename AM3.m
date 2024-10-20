function [LHS,RHS] = AM3(timestep,M,B,C,sol)
arguments
    timestep (1,1) double
    M (:,:) double
    B (:,:) double
    C (:,3) double
    sol (:,2) double
end
%AM3 Adams-Moulton-method of 3rd order for solving linear inhomogeneous
%differential equations.
% The equation should have the following layout:
%
%     ∂ϕ
%   M*-- = B(t)*ϕ + C(t)
%     ∂t
%
% Input:
%   timestep: Integration interval
%   M:        Factor in front of dϕ/dt
%                M ∈ m×m
%   B:        Factor in front of ϕ
%                B ∈ m×3m   Layout: [B(t_n+1), B(t_n), B(t_n-1)]
%   C:        Inhomogeneous part
%                C ∈ m×3    Layout: [C(t_n+1), C(t_n), C(t_n-1)]
%   sol:      solution of the timestep(s) before
%                sol ∈ m×2  Layout: [ϕ(t_n), ϕ(t_n-1)]
%
% Output:
%   LHS: Left-hand-side of the resulting equation system
%   RHS: Right-hand-side of the resulting equation system
%
% Exercise 6
%
% © 2024, Andreas Steger

% get width of matrices
n = width(M); % must be width

% check the dimensions
if (n ~= width(M))
    error("M must be a square matrix (is "+n+"×"+width(M)+")")
end
if (width(B) ~= 3*height(B))
    error("B must be composed of three square matrices (all three currently " + ...
        "interpreted as "+height(B)+"×"+width(B)/3+")")
end
if (n ~= height(B))
    error("M and B must be of same height (is "+height(B)+", expected: "+n+")")
end
if(height(C) ~= n)
    error("C must be composed of row vectors with height of M")
end
if(height(sol) ~= n)
    error("sol must be composed of row vectors with height of M")
end

% do the calculation
LHS = M - (5*timestep)/12 * B(:,1:n);
RHS = M*sol(:,1) + timestep/12*(...
    5*C(:,1)+...
    8*(B(:,n+1:2*n) * sol(:,1) + C(:,2))-...
    (B(:,2*n+1:3*n) * sol(:,2) + C(:,3)));
end