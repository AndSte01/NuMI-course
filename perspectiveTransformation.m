function x = perspectiveTransformation(nodes_source, nodes_destination, point)
arguments
    nodes_source (4,2) double
    nodes_destination (4,2) double
    point (2,1) double
end
%PERSPECTIVETRANSFORMATION transforms a point from one 2D-Quadrangle
%area to another. It uses perspective transformation to do so.
%   Note: this function has severe rounding issues if one of the forms 
%   differs to much from a rectangular shape.
%
%   Node-lists have the following from:
%       [x, y; ...; x, y]
%
%                  _-3                                   ____----3'
%               _-‾   \                      ____----‾‾‾‾       /
%            _-‾       \           4'----‾‾‾‾                 /
%         4-‾            \          \                        /
%        /     'from'     \   -->    \        'to'         /
%       /              __--2          \                  /
%      /         __--‾‾                \                /
%     /    __--‾‾                      1'---___       /
%    1--‾‾                                     ‾‾‾---2'
%
% Input:
%   nodes_source:      the nodes of the 2D-Quadrangle to transform FROM
%   nodes_destination: the nodes of the 2D-Quadrangle to transform to
%   point:             the point in the source quadrangle that shall be
%                      transformed [x; y]
%
% Output:
%   x: the transformed point
%
% Exercise 5
%
% © 2024, Andreas Steger

% numerical solutions derived from the formulas given at using Matlab
% symbolic calculation
% https://math.stackexchange.com/questions/296794/finding-the-transform-matrix-from-4-projected-points-with-javascript/339033#339033
A_inv_num = @(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4)...
    [(y_2 - y_3)/(x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3), -(x_2 - x_3)/(x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3), (x_2*y_3 - x_3*y_2)/(x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3);...
    (y_1 - y_3)/(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3), -(x_1 - x_3)/(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3), (x_1*y_3 - x_3*y_1)/(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3);...
    (y_1 - y_2)/(x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2), -(x_1 - x_2)/(x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2), (x_1*y_2 - x_2*y_1)/(x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2)];
B_num = @(x_1, y_1, x_2, y_2, x_3, y_3, x_4, y_4)...
    [(x_1*(x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), -(x_2*(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), (x_3*(x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2);...
    (y_1*(x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), -(y_2*(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), (y_3*(x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2))/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2);...
    (x_2*y_3 - x_3*y_2 - x_2*y_4 + x_4*y_2 + x_3*y_4 - x_4*y_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), -(x_1*y_3 - x_3*y_1 - x_1*y_4 + x_4*y_1 + x_3*y_4 - x_4*y_3)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2), (x_1*y_2 - x_2*y_1 - x_1*y_4 + x_4*y_1 + x_2*y_4 - x_4*y_2)/(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2)];

% do the transformation
x_tmp_num = B_num(nodes_destination(1, 1), nodes_destination(1, 2), nodes_destination(2, 1), nodes_destination(2, 2), ...
    nodes_destination(3, 1), nodes_destination(3, 2), nodes_destination(4, 1), nodes_destination(4, 2))...
    *A_inv_num(nodes_source(1, 1), nodes_source(1, 2), nodes_source(2, 1), nodes_source(2, 2), ...
    nodes_source(3, 1), nodes_source(3, 2), nodes_source(4, 1), nodes_source(4, 2))...
    *[point; 1];

% calculate the values of the new point
x = [x_tmp_num(1) x_tmp_num(2)]'/x_tmp_num(3);

end