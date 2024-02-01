function x = getxPosAna(nodes, xi, eta)
arguments
    nodes (4,2) double
    xi (1,1) double
    eta (1,1) double
end
%GETXPOS Maps xi and eta to the 2D-Space provided by given nodes
% This function uses the analytical solution
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
    x(:,i) = nodes(1,:) + (nodes(2,:)-nodes(1,:)) * (xi(i)+1)/2 +...
    ((nodes(4,:)-nodes(1,:)) + ((nodes(3,:)-nodes(4,:))-(nodes(2,:)-nodes(1,:))) * (xi(i)+1)/2)...
    * (eta(i)+1)/2;
end
end