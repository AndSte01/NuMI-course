function [points] = extractElementPointsFromList(x, ele)
arguments
    x (:,2) double
    ele (1,:) double
end
%EXTRACTELEMENTPOINTSFROMLIST Extracts the points given by the indices
%ele in x
%
% Input:
%   x:   The points
%            [x, y; ...; x, y]
%   ele: The indices in x of the points (stored at that indices in x)
%
% Output:
%   points: The extracted points
%               [x, y; ...; x, y]
%
% Helper for Exercise 7, 8, 10
%
% © 2024, Andreas Steger

%find out how many points we want to get
n_points = length(ele);

% empty result vector
points = zeros([n_points 2]);

% now extract the points
for i=1:n_points
    points(i,:) = x(ele(i),:);
end
end