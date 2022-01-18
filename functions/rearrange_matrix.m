function output = rearrange_matrix(data,row,dim)
%REARRANGE_MATRIX rearranges a matrix such that row becomes the first row,
%while maintaining order.
%   data: matrix you want rearranged
%   row: the row at which you want to become the first
if dim == 3
    n = size(data);
    mat1 = data(1:row-1,:,:);
    mat2 = data(row:n,:,:);
    output = [mat2; mat1]; 
elseif dim == 4
    n = size(data);
    mat1 = data(1:row-1,:,:,:);
    mat2 = data(row:n,:,:,:);
    output = [mat2; mat1]; 
else
    n = size(data);
    mat1 = data(1:row-1,:);
    mat2 = data(row:n,:);
    output = [mat2; mat1];
end
end

