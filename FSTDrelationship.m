% Compare the FSD of different ice thicknesses


%% Read in the data.
clear
close all
addpath functions
addpath packages/joyPlot

filedir = "cases/init/history/iceh.2005-01-01.nc";
           
filecase = "iceh.2005-01-01";

ncdisp(filedir)

% Parameters
fstd_switch = 0;
coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
grid = 'gx1';
[lat,lon,row] = grid_read(grid);

coords = [-60,20];
data = data_format(filedir,'aicen_d',row,lat,lon,3); %aicen
[lat_out,lon_out] = lat_lon_finder(coords(1),coords(2),lat,lon);
data(lon_out,lat_out,2)

data = data_format(filedir,'afsdn_d',row,lat,lon,4); %aicen

NFSD = ncread(filedir,'NFSD');
NCAT = ncread(filedir,'NCAT');
% Cell 
%figure(1)
t1=tiledlayout(1,5);

for i = 1:5
    data_temp = zeros(1,24);
    for j = 1:24
        norm_FSD(i,j) = data(lon_out,lat_out,j,i);
    end
    %norm_FSD = frac_pdf(data_temp);
    x_axis = 1:length(norm_FSD);
    %if sum(norm_FSD(i,:)) > eps
%         nexttile(t1)
%         area(x_axis,norm_FSD(i,:))
%         title(strcat(sprintf("FSD from ITD cat %d, ",i), filecase))
%         ylabel("probability")
%         xlabel("FSD bins")
%         ylim([0,0.005])
    %end
end
norm_FSD1 = norm_FSD(1,:);
figure(1)
area(1:24,norm_FSD1)
xlabel("FSD categories")
ylabel("ITD categories")
title("Cell 1st cat")
norm_FSD1 = norm_FSD(1,:);
norm_FSDrest = norm_FSD(2:end,:);
figure(2)
joyPlot(norm_FSDrest',1:24,0.005,'FaceColor',[0 0.4470 0.7410],'StrokeColor','w')
xlabel("FSD categories")
ylabel("ITD categories")
title("Cell rest")

%%
figure(3)
t2=tiledlayout(3,8);

for i = 1:24
    data_temp = zeros(1,5);
    for j = 1:5
        norm_ITD(i,j) = data(lon_out,lat_out,i,j);
    end
    x_axis = 1:length(norm_ITD(i,:));
%     if sum(norm_ITD(i,:)) > eps
%         area(x_axis,norm_ITD(i,:))
%         title(strcat(sprintf("ITD from FSD cat %d, ",i), filecase))
%         ylabel("probability")
%         xlabel("ITD bins")
%     end
end
figure(2)
joyPlot(norm_ITD',1:5,0.005,'FaceColor',[0 0.4470 0.7410],'StrokeColor','w')
xlabel("ITD categories")
ylabel("FSD categories")
title("Cell")

data = data_format(filedir,'afsd_d',row,lat,lon,3);
for j = 1:24
    norm_AFSD(j) = data(lon_out,lat_out,j);
end

figure(3)
area(1:24,norm_AFSD)
xlabel("FSD categories")
ylabel("ITD categories")
title("Cell")



% Global
data = data_format(filedir,'afsdn_d',row,lat,lon,4);
% figure(4)
% t3=tiledlayout(3,8);
% norm_ITD = frac_pdf(data,1);
% x_axis = 1:5;
% for i = 1:24
%     if sum(norm_ITD(i,:)) > eps
%         nexttile(t3)
%         area(x_axis,norm_ITD(i,:))
%         title(sprintf("Global ITD from FSD cat %d, ",i))
%         ylabel("probability")
%         xlabel("ITD bins")
%     end
% end
figure(4)
joyPlot(norm_ITD',1:5,0.005,'FaceColor',[0 0.4470 0.7410],'StrokeColor','w')
xlabel("ITD categories")
ylabel("FSD categories")
title("Global")


figure(5)
% t4=tiledlayout(3,8);
% 
% norm_FSD = frac_pdf(data,2);
% x_axis = 1:24;
% for i = 1:5
%     if sum(norm_FSD(i,:)) > eps
%         nexttile(t4)
%         area(x_axis,norm_FSD(i,:))
%         title(sprintf("Global FSD from ITD cat %d, ",i))
%         ylabel("probability")
%         xlabel("FSD bins")
%         %count = count + 1;
%     end
% end
joyPlot(norm_FSD',1:24,0.005,'FaceColor',[0 0.4470 0.7410],'StrokeColor','w')
xlabel("FSD categories")
ylabel("ITD categories")
title("Global")

