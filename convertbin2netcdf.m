clear all;close all;clc;
f = fopen('grid_gx1.bin');

z1 = [];
for i = 1:320
    for j = 1:384
        a(i,j) = fread(f,1,'float');
    end
    z(:,:) = a(i,:);
    z1 = [z1;z];
    clear z
end
