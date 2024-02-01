function [x,i] = solveG(A,b,x0,rtol,itermax)
arguments
    A (:,:) double
    b (:,1) double
    x0 (:,1) double
    rtol (1,1) double
    itermax (1,1) int64
end
%SOLVEG Solves the linear system A*x = b for x using the iterative gradient
%approach
%
%   A*x = b
%
% Input:
%   A:       The matrix A
%                A ∈ n×n
%   b:       The vector b
%                b ∈ n×1 [b_1; ...; b_n] (column vector)
%   x0:      Starting position of the iterator
%                x0 ∈ n×1 [x_1; ...; x_n] (column vector)
%   rtol:    minimal error that is acceptable (iteration stops)
%   itermax: maximal number of iterations before the algorithm stops
%
% Output:
%   x: The solution
%          x ∈ n×1
%   i: The number of iterations
%
% © 2024, Andreas Steger

% test inputs for validity
size = height(A);
if(height(A) ~= width(A))
    error("A must be a square matrix (is "+size+"×"+width(A)+")")
end
if(height(b)~=size)
    error("b must be vector of same height as A (is "+height(B)+" and "+size+")")
end
if(height(x0)~=size)
    error("x0 must be vector of same height as A (is "+height(B)+" and "+size+")")
end

% set initial values
i    = 0;       % counting variable
x    = x0;      % start value (given by user)
r    = b - A*x; % initial direction
rcur = norm(r); % current error (make sure the initial value is bigger)

% do the iterations
while i < itermax && rcur > rtol
    % increase iteration counter
    i = i+1;
    
    % calculate new x
    v = A*r;
    alpha = (r'*r)/(r'*v);
    x = x + alpha*r;
    r = r - alpha*v;
    
    % calculate the new error
    rcur = norm(r);
end

end