function x = getxPos(nodes, xi, eta)
arguments
    nodes (4,2) double
    xi (1,1) double
    eta (1,1) double
end
%GETXPOS Maps xi and eta to the 2D-Space provided by given nodes
% This function uses a numerical solution
%
%                                        _-3'
%  [-1,1]        [1,1]                _-‾   \
%     4 ---------- 3               _-‾       \
%     |            |            4'‾            \
% eta | Omega_ref  |   -->     /     Omega      \
%     |            |          /              __--2'
%     |            |         /         __--‾‾
%     1 -----------2        /    __--‾‾
%  [-1,-1]      [1,-1]     1'--‾‾
%
% Input:
%   nodes: Nodes of the 2D-Space (oder see diagram above)
%              layout: [x_1, y_1; ...; x_4, y_4]
%   xi:    fist value of coordinate in reference system
%   eta:   second value of coordinate in reference system
%
% Output:
%   x: Coordinates in original system
%
% Exercise 5
%
% requires linquadref
%
% © 2024, Andreas Steger

% some error checking
if (length(xi) ~= length(eta))
    error("lengths of xi and eta must match");
end

% number of points to calculate
n = length(xi);

% calculate x
x = zeros([2 n]);

for i = 1:n
    x(:,i) = linquadref(xi(i), eta(i))' * nodes;
end
end