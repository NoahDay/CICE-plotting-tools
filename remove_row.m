function output = remove_row(data, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
mat1 = data(1:n-1,:);
mat2 = data(n+1:size(data),:);
output = [mat1;mat2]; 
end

