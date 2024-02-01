function df = finDifSel(method, h, grid_points)
arguments
    method (1,1) string {mustBeMember(method,["2ptLeft","2ptRight","3ptRight","3ptCenter","5ptCenter"])}
    h (1,1) double
    grid_points (:,1) double
end
%finDifSel calculates the first derivative using different methods of finite
%differences
%
% The following methods are supported:
%   2ptLeft     [x0 − h, x0]
%   2ptRight    [x0, x0 + h]
%   3ptRight    [x0, x0 + h, x0 + 2h]
%   3ptCenter   [x0 − h, x0 + h]
%   5ptCenter   [x0 − 2h, x0 − h, x0 + h, x0 + 2h]
%
% Input:
%   method:      string selecting one of the above methods
%   h:           step width of the gird
%   grid_points: the points of the function at the above documented points
%
% Output:
%   df: the calculated value of the derivative
%
% Helper for Exercise 4.1
%
% © 2024, Andreas Steger

switch method
    case {"2ptLeft", "2ptRight"}
        df = (grid_points(2) - grid_points(1)) / h;

    case "3ptRight"
        df = (-3*grid_points(1) + 4*grid_points(2) - grid_points(3))/(2*h);

    case "3ptCenter"
        df = (grid_points(2) - grid_points(1))/(2*h);

    case "5ptCenter"
        df = (grid_points(1) - 8*grid_points(2) + 8*grid_points(3) - grid_points(4))/(12*h);

    otherwise
        error("The selected method '" + method +"' is not supported/implemented");
end
end