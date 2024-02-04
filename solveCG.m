function [x,i] = solveCG(A,b,x0,rtol,itermax)
arguments
    A (:,:) double
    b (:,1) double
    x0 (:,1) double
    rtol (1,1) double
    itermax (1,1) int64
end
%SOLVEG Solves the linear system A*x = b for x using the iterative
%conjugate gradient approach
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
    error("b must be vector of same height as A (is "+height(B)+", expected: "+size+")")
end
if(height(x0)~=size)
    error("x0 must be vector of same height as A (is "+height(B)+", expected: "+size+")")
end

% set initial values
i    = 0;       % counting variable
x    = x0;      % start value (given by user)
r    = b - A*x; % initial direction
p    = r;
rcur = norm(r); % current error (make sure the initial value is bigger)

% do the iterations
while i < itermax && rcur > rtol
    % increase iteration counter
    i = i+1;
    
    % calculate new x
    v = A*p;
    alpha = (r'*r)/(p'*v);
    x = x + alpha*p;
    r_n = r - alpha*v;
    beta = (r_n'*r_n)/(r'*r);
    p = r_n + beta*p;
    r = r_n; % 'forget' the old r

    % calculate the new error
    rcur = norm(r);
end

end