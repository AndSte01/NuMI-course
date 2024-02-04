function x = solveGauss(A,b)
arguments
    A (:,:) double
    b (:,1) double
end
%SOLVEGAUSS Solves a linear System A*x = b for x using the Gauss approach
%
%   A*x = b
%
% Input:
%   A: The matrix A
%          A ∈ n×n
%   b: The vector b
%          b ∈ n×1 [b_1; ...; b_n] (column vector)
%   
% Output:
%   x: The solution
%          x ∈ n×1
%
% © 2024, Andreas Steger

% get the number of of rows/columns (A is quadratic)
size = height(A);

% test inputs for validity
if(height(A) ~= width(A))
    error("A must be a square matrix (is "+size+"×"+width(A)+")")
end
if(height(b)~=size)
    error("b must be vector of same height as A (is "+height(B)+", expected: "+size+")")
end

% iterate over columns
for c=1:(size-1)
    % iterate over rows
    for r=(c+1):size
        fac = A(r,c)/A(c,c); % Prefactor
        A(r,:) = A(r,:)-fac*A(c,:);
        b(r) = b(r)-fac*b(c);
    end
end

% empty vector for x
x = zeros([size 1]);

% backward substitution
for r=-size:-1 % why MatLab why? (count down not up)
    % note: x = zeros -> only important values are not null
    x(-r) = (b(-r) - A(-r,:)*x)/A(-r,-r);
end

% for possible implementation of inversion
% A(3,:) = A(2,:)-A(2,3)/A(3,3)*A(3,:);

end