function done = patch_data(data,row,copies)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n = size(data);
mat1 = data(1:row-1,:);
mat2 = data(row:n,:);
insert_vec = (data(row-1,:) + data(row,:))/2;

done = [mat1; repmat(insert_vec,copies,1); mat2];
end

