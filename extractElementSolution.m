function [T_ele] = extractElementSolution(T, ele)
arguments
    T (1,:) double
    ele (1,:) double
end
%EXTRACTELEMENTSOLUTION Extracts the solution only for the points in the
%element.
%
% Input:
%   T:   The global solution
%   ele: The indices in x of the points (stored at that indices in x)
%
% Output:
%   T_ele: The extracted solution
%               [T; ...; T]
%
% Helper for Exercise 8, 10
%
% © 2024, Andreas Steger

%find out how many points we want to get
n_points = length(ele);

% empty result vector
T_ele = zeros([n_points 1]);

% now extract the points
for i=1:n_points
    T_ele(i) = T(ele(i));
end

end