function [done, n] = patch_mat(matrix)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[n, m] = size(matrix);

idx = matrix(:,1) == 0;

row = find((1:n).*idx');

mat1 = matrix(1:row-1,:);
mat2 = matrix(row:n,:);
insert_vec = (matrix(row-1,:) + matrix(row,:))/2;

done = [mat1; insert_vec; mat2];
end

