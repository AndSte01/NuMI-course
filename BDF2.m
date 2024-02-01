function [LHS,RHS] = BDF2(timestep,M,B,C,sol)
arguments
    timestep (1,1) double
    M (:,:) double
    B (:,:) double
    C (:,1) double
    sol (:,2) double
end
%BDF2 BDF2-method for solving linear inhomogeneous differential equations.
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
%                B ∈ m×m    Layout: B(t_n+1)
%   C:        Inhomogeneous part
%                C ∈ m×1    Layout: C(t_n+1)
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

% check the dimensions
n = height(M);
if (n ~= width(M))
    error("M must be a square matrix (is "+n+"×"+width(M)+")")
end
if (n ~= height(B))
    error("M and B must be of same size ( M ∈ "+n+"×"+width(M)+" and B ∈ "+ ...
        height(B)+"×"+width(B)+")")
end
if (height(C) ~= n)
    error("C must be composed of row vectors with the same height as M")
end
if(height(sol) ~= n)
    error("sol must be row vector with height of M")
end

% do the calculation
LHS = 3/2*M - timestep*B;
RHS = 2*M*sol(:,1) - 1/2*M*sol(:,2) + timestep*C;
end