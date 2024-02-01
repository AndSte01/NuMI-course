function [sysmat,rhs] = assemble(elemat,elevec,sysmat,rhs,ele)
arguments
    elemat (:,:) double
    elevec (:,1) double
    sysmat (:,:) double
    rhs (:,1) double 
    ele (:,1)
end
%ASSEMBLE Assemble an element matrix into the system matrix
%
% Input:
%   elemat: The element Matrix
%   elevec: The element vector
%   sysmat: The system matrix
%   rhs:    The system vector
%   ele:    global node indices
%               [n_1; n_2; ...; n_N]
%
% Output:
%   sysmat: The updated system matrix
%   rhs:    The updated system vector
%
% Exercise 7
%
% © 2024, Andreas Steger

% get the number of nodes
n_N = length(ele);

% some basic error checking
if n_N ~= 4
    error("number of nodes must be equal to 4");
end

% iterate over the element matrix and assemble accordingly
for i=1:n_N
    for j=1:n_N
        sysmat(ele(i),ele(j)) = sysmat(ele(i),ele(j)) + elemat(i,j);
    end

    % while we are here assemble the system vector as well
    rhs(ele(i)) = rhs(ele(i)) + elevec(i);
end

end