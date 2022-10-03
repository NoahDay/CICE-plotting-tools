% 3 cases:
%   1. With full WIM
%   2. With WIM off
%   3. With Waves off


% Setup
clear all
clc
close all

%%
addpath functions
addpath /Users/noahday/GitHub/CICE-analyser/processing

historydir.waves = '/Volumes/NoahDay5TB/WIMonAlessandroRun/history/2019/';
historydir.sheet = '/Volumes/NoahDay5TB/WIMoffSheetIceAlessandroRun/history/';
historydir.pancake = '/Volumes/NoahDay5TB/WIMoffPancakeAlessandroRun/history/';


a.waves = dir([historydir.waves '/*.nc']);
n_files = numel(a.waves);

for i = 1:n_files
   filenames.waves(i,:) = strcat(historydir.waves,a.waves(i).name);
   dirdates.waves(i,:) = a.waves(i).name(6:end-3);
end

a.sheet = dir([historydir.sheet '/*.nc']);
n_files = numel(a.sheet);
for i = 1:n_files
   filenames.sheet(i,:) = strcat(historydir.sheet,a.sheet(i).name);
   dirdates.sheet(i,:) = a.sheet(i).name(6:end-3);
end

a.pancake = dir([historydir.pancake '/*.nc']);
n_files = numel(a.pancake);
for i = 1:n_files
   filenames.pancake(i,:) = strcat(historydir.pancake,a.pancake(i).name);
   dirdates.pancake(i,:) = a.pancake(i).name(6:end-3);
end

%% Select the relevant variables
% Columns are data
sector = "SH";
close all
var_list = {'aice','hi','iage','fsdrad','frazil','frzmlt','congel','sice','strint','dafsd_latg','dafsd_latm','dafsd_newi','dafsd_weld','dafsd_wave','full_fstd'};
[X_temp, row_idx]= read_data_vec(filenames.sheet,sector,var_list); % [var_list, lon, lat]
%[~, ~,dep] = size(X_temp);

%clear X_temp
label_vec = variable_dict(var_list);
size(X_temp)
%
data.Xunstandard = X_temp;
data.row_idx = row_idx;
save_filename = strcat('X_sheet.mat');
save(save_filename,'data','-v7.3');
clear data
%%
clear X row_idx
load('X_waves.mat')
X.waves = data.Xunstandard;
row_idx.waves = data.row_idx;
clear data

load('X_sheet.mat')
X.sheet = data.Xunstandard;
row_idx.sheet = data.row_idx;
clear data

load('X_pancake.mat')
X.pancake = data.Xunstandard;
row_idx.pancake = data.row_idx;
clear data

%% Clean the data
% NEED TO DO THIS ALL TOGETHER SO ITS ALL STANDARDISED CORRECTLY

[Xnan.waves,row_idx.waves] = clearNaN(X.waves);
disp('Waves done!')
[Xnan.sheet, row_idx.sheet] = clearNaN(X.sheet);
disp('Sheet done!')
[Xnan.pancake, row_idx.pancake] = clearNaN(X.pancake);
disp('Pancake done!')

dimension.waves = size(Xnan.waves);
dimension.sheet = size(Xnan.sheet);
dimension.pancake = size(Xnan.pancake);
Xnan.all = [Xnan.waves;Xnan.sheet;Xnan.pancake];
%%
cleaned_data.X_no_nan = Xnan.all;
cleaned_data.row_idx = row_idx;
cleaned_data.dimension = dimension;
save('cleanded_kmeans_data.mat','cleaned_data');
%%
%load('kmeans_3cases_4_classes.mat')
%% Standardise the data
X_standard_all = Xnan.all;
[len,wid] = size(X_standard_all);
% DONT STANDARDISE THE LAT AND LON
for j = 1:wid
    if min(Xnan.all(:,j)) < 0
        % If negatives, then recentre so its all positive
        max_X(j) = max(X_standard_all(:,j));
        min_X(j) = min(X_standard_all(:,j));
        X_standard_all(:,j) = (X_standard_all(:,j) - min_X(j))/(max_X(j) - min_X(j));
    end
    %Step 1: Log transformation as almost all of the data is highly skewed
    X_standard_all(:,j) = log(X_standard_all(:,j)+1);
    %Step 2: Standardization
    % Calculate mean
    mean_X(j) = mean(X_standard_all(:,j));
     % Calculate standard deviation
    std_X(j) = std(X_standard_all(:,j));
    
    X_standard_all(:,j) = X_standard_all(:,j)/std_X(j);

    % Standardise the data
   max_X(j) = max(X_standard_all(:,j));
   min_X(j) = min(X_standard_all(:,j));
   X_standard_all(:,j) = (X_standard_all(:,j) - min_X(j))/(max_X(j) - min_X(j));
end
%clear X_temp
%%
eva = evalclusters(X_standard_all(:,1:end-2),'kmeans','CalinskiHarabasz','KList',1:10);
temp = eva.CriterionValues;
f = figure;
plot(eva)
%% k-means clustering


X = X_standard_all(:,1:end-2);
num_clusters = 4;
[idx,C] = kmeans(X,num_clusters,'MaxIter',300);
kmeans_cluster.idx = idx;
kmeans_cluster.row_idx = row_idx;
kmeans_cluster.label_vec = label_vec;
kmeans_cluster.C = C;
kmeans_cluster.X_standard_all = X_standard_all;
kmeans_cluster.num_clusters = num_clusters;
save_filename = strcat('kmeans_3cases_',sprintf('%g',num_clusters),'_classes.mat');
save(save_filename,'kmeans_cluster');

%% Average stats
X_temp = [];
[~,~,dep] = size(Xnan.all);
for i = 1:dep
    X_temp = [X_temp; Xnan.all(:,:,i)];
end

X_temp = [X_temp(:,1:end-2), idx];

for i = 1:num_clusters
    index_class_data = X_temp(:,end) == i;
    class_data = X_temp(index_class_data,1:end-1);
    average_stats(i,:) = mean(class_data);
end

close all
conFigure(11)
f = figure('Position',[100, 100, 100, 100]);
%set(gcf,'Visible', 'on')
for i = 1:13
  subplot(7,2,i)
  bar(average_stats(:,i))
  xticks(1:length(average_stats(:,i)))
  %set(gca, 'XTickLabel',label_vec)
  title(label_vec{i})
end
exportgraphics(f,'statcomp.pdf','ContentType','vector')

%% Distributions of each class

close all
conFigure(30)
f = figure('Position',[100, 100, 100, 100]);
%set(gcf,'Visible', 'on')
for i = 1:num_clusters
index_class_data = X_temp(:,end) == i;
class_data = X_temp(index_class_data,3);
  subplot(2,2,i)
  hist(class_data)
  xticks(1:length(average_stats(:,i)))
  %set(gca, 'XTickLabel',label_vec)
  title(sprintf('Class %g',i))
end
exportgraphics(f,'age_dist.pdf','ContentType','vector')

%%
index_lat = X_new(:,end-1) == X_new(1,end-1);
index_lon = X_new(:,end) == X_new(1,end);
index_both = index_lat.*index_lon;
sum(index_both)

index_lat = X_new(:,end-1) == X_new(10,end-1);
index_lon = X_new(:,end) == X_new(5,end);
index_both = index_lat.*index_lon;
sum(index_both)
%%
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;

addpath functions
[lat,lon,~,ulat,ulon] = grid_read('om2');
clear X_map MIZ_width
% Make worldmap with colour matching
uarea = data_format(filenames.waves(1,:),"uarea");

[num_files,~] = size(filenames.waves);
sector = "EA";
coords = sector_coords(sector);
font_size = 5;
plot_type = "kmeans";
plotting = "on";
conFigure(11)
X_standard_all(:,end-1:end) = Xnan.all(:,end-1:end);
X_new = X_standard_all(1:dimension.waves(1),1:dimension.waves(2));
[len,wid] = size(lat);

file_number_vec = 1:365;
for file_number = file_number_vec %num_files-1
    row_file = [0,cumsum(row_idx.waves)];
    if plot_type == "pc1"
        file_idx = score(row_file(file_number)+1:row_file(file_number+1),1);
        file_idx = file_idx + abs(min(file_idx));
        file_idx = file_idx./max(file_idx);
    elseif plot_type == "kmeans"
        index_lat = X_new(:,end-1) == X_new(10,end-1);
        index_lon = X_new(:,end) == X_new(10,end);
        index_both = index_lat.*index_lon;clc
        diff_index = diff(find(index_both)');
        row_file = [0,cumsum(diff_index)];
        %row_file = find(index_both);%[0,cumsum(row_idx)];
        row_vec = row_file(file_number)+1:row_file(file_number+1);
        file_idx = idx(row_vec);
    end
    %[aice, sector_mask] = data_format_sector(filenames.waves(file_number,:),'aice',sector);
    %[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    
        %ice_mask =  aice > 0.15;
        %k_means(~ice_mask) = NaN;
    if plotting == "on"
        f = figure;
        set(gcf,'Visible', 'off')
        w = worldmap('world');
            %axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            %setm(w, 'Origin', [-90 0 0]);
            %setm(w, 'maplatlimit', [-90,-50]);
            
            %setm(w, 'maplonlimit', [coords(1,2),coords(3,2)]);

            axesm eqdcylin; 
            setm(w, 'Origin', [0 28 0]); 
            setm(w, 'maplatlimit', [-75,-50]); setm(w, 'maplonlimit', [1,150]); 

            setm(w, 'meridianlabel', 'off')
            setm(w, 'parallellabel', 'off') 
            setm(w, 'mlabellocation', 60);
            setm(w, 'plabellocation', 10);
            setm(w, 'mlabelparallel', -45);
            setm(w, 'grid', 'off');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,k_means)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
            colorbar; cmocean('deep');
            title(dirdates.waves(file_number,:),'Color','black','FontSize',font_size+5)
            %plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
            
            if plot_type == "pc1"
                caxis([0,1])
            elseif plot_type == "kmeans"
                cb = colorbar; cmocean('deep',num_clusters)
                %cb.TickLabels = region_label; 
                cb.Ticks = 1:num_clusters;
                cb.Location = 'southoutside';%cb.Location = 'northoutside';
                caxis([1,num_clusters])
            end
    end
    % Calculate the width of the MIZ
    % CHECK WHAT number MIZ is !!!!!!!!
    [MIZ_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames.waves(file_number,:)),sector,k_means,1);
    

    % Calculate the perimeter and area of the zones
    for rgn_number = 1:num_clusters
        region_data = k_means == rgn_number;
        [len, wid] = size(k_means);
        perimeter_region(rgn_number,file_number) = find_perimeter(region_data,len,wid);

        % Area
        area_region(rgn_number, file_number) = sum(sum(uarea.*region_data));
    end
    if plotting == "on"
        if file_number == file_number_vec(1)
            gif('kmeans_EA_SA_noMinSIC_stndzd_4class.gif','DelayTime',0.5,'resolution',500,'overwrite',true)
        else
            gif
        end
    end
    clear aice ice_mask k_means file_idx X_map
end

%%
region_label = {'1','2','3','4'};
date_label = datetime(dirdates.waves);
date_label_vec = date_label(file_number_vec);
conFigure(11)
figure
subplot(2,3,1)
plot(date_label_vec,[perimeter_region; sum(perimeter_region)],'LineWidth',5)
ylabel('Perimeter [block edges]')
%legend(region_label,'Location','southeast')

subplot(2,3,2)
plot(date_label_vec,[area_region.*1e-6.*1e-6; sum(area_region.*1e-6.*1e-6)],'LineWidth',5)
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')

% Create a tile on the right column to get its position
legend_label = region_label;
legend_label{end+1} = 'Sum of regions';
ax = subplot(2,3,3,'Visible','off');
axPos = ax.Position;
delete(ax)
% Construct a Legend with the data from the sub-plots
hL = legend(legend_label);
% Move the legend to the position of the extra axes
hL.Position(1:2) = axPos(1:2);

subplot(2,3,4)
perim_per_area = perimeter_region./(area_region.*1e-6.*1e-6);
plot(date_label_vec,perim_per_area,'LineWidth',5)
ylabel('Perimeter per Area [block per million km$^2$]')
%text(1,1,sprintf("$sigma =$ %g",std(perimeter_region'./(area_region./1000)')))
% ta = annotation('textarrow', [0.76 0.83], [0.71 0.71]);
% ta.String = sprintf("$\sigma =$ %g",std(perimeter_region(4,:)'./(area_region(4,:)./1000)'));
% ta.Interpreter = "latex";
%legend(region_label,'Location','southeast')



diff_perim = perim_per_area(:,2:end)' - perim_per_area(:,1:end-1)';
subplot(2,3,5)
bar(std(diff_perim))
xticklabels(region_label)
ylabel('S.D. $\Delta$ perimeter per area')










%% FUNCTIONS
function [X_total, idx] = read_data_vec(filenames,sector,var_list)
    [len,~] = size(filenames);
    if len == 1
        filename = filenames;
    else
        filename = filenames(1,:);
    end
    
    SIC_threshold = 0.01;
    [aice, sector_mask] = data_format_sector(filename,"aice",sector);
    [lat,lon,~,ulat,ulon] = grid_read('om2');
    lat_vec = reshape(lat(sector_mask),numel(lat(sector_mask)),1);
    lon_vec = reshape(lon(sector_mask),numel(lon(sector_mask)),1);

    NFSD = ncread(filename,"NFSD");
    NCAT = ncread(filename,"NCAT");
    [floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
    
    X_total = [];
    for j = 1:len
        filename = filenames(j,:);
        [aice, sector_mask] = data_format_sector(filename,"aice",sector);
        ice_mask = aice > SIC_threshold;
        X_temp = [];
        for i = 1:numel(var_list)
            if var_list{i} == "pancake"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp2(ii,jj) =  temp(ii,jj,1,1)./aice(ii,jj).*floe_binwidth(1).*100;
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "full_fsd"
                temp = data_format_sector(filename,"afsd",sector);
                [len, wid] = size(aice);
                for kk = 1:length(NFSD)
                    for ii = 1:len
                        for jj = 1:wid
                            if aice(ii,jj) > SIC_threshold
                                temp2(ii,jj) =  temp(ii,jj,kk)./aice(ii,jj).*floe_binwidth(kk).*100;
                            else
                                temp2(ii,jj) =  NaN;
                            end
                        end
                    end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
                end
            elseif var_list{i} == "full_fstd"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for nc = 1:length(NCAT)
                for kk = 1:length(NFSD)
                    for ii = 1:len
                        for jj = 1:wid
                            if aice(ii,jj) > SIC_threshold
                                temp2(ii,jj) =  temp(ii,jj,kk,nc)./aice(ii,jj).*floe_binwidth(kk).*100;
                            else
                                temp2(ii,jj) =  NaN;
                            end
                        end
                    end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
                end
                end
            elseif var_list{i} == "thick_pan"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            for ll = 2:5
                                temp2(ii,jj) =  temp2(ii,jj) + temp(ii,jj,1,ll)./aice(ii,jj).*floe_binwidth(1).*100;
                            end
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                %idx = isnan(temp2);
                %temp2(idx) = 0;
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "large"
                temp = data_format_sector(filename,"afsdn",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp2(ii,jj) =  temp(ii,jj,numel(NFSD),1)./aice(ii,jj).*floe_binwidth(numel(NFSD)).*100;
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];

            elseif var_list{i} == "dafsd_newi"
                temp = data_format_sector(filename,"dafsd_newi",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_weld"
                temp = data_format_sector(filename,"dafsd_weld",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_wave"
                temp = data_format_sector(filename,"dafsd_wave",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "dafsd_latm"
                temp = data_format_sector(filename,"dafsd_latm",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
             elseif var_list{i} == "dafsd_latg"
                temp = data_format_sector(filename,"dafsd_latg",sector);
                [len, wid] = size(aice);
                for ii = 1:len
                    for jj = 1:wid
                        if aice(ii,jj) > SIC_threshold
                            temp_fsd(:) = temp(ii,jj,:);
                            temp2(ii,jj) =  sum((temp_fsd./aice(ii,jj)).*NFSD');
                        else
                            temp2(ii,jj) =  NaN;
                        end
                    end
                end
                temp_vec = reshape(temp2(sector_mask),numel(temp2(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strair"
                tempx = data_format_sector(filename,"strairx",sector);
                tempy = data_format_sector(filename,"strairy",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strocn"
                tempx = data_format_sector(filename,"strocnx",sector);
                tempy = data_format_sector(filename,"strocny",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strint"
                tempx = data_format_sector(filename,"strintx",sector);
                tempy = data_format_sector(filename,"strinty",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "strcor"
                tempx = data_format_sector(filename,"strcorx",sector);
                tempy = data_format_sector(filename,"strcory",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];
            elseif var_list{i} == "vel"
                tempx = data_format_sector(filename,"uvel",sector);
                tempy = data_format_sector(filename,"vvel",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];  
            elseif var_list{i} == "atmspd"
                tempx = data_format_sector(filename,"uatm",sector);
                tempy = data_format_sector(filename,"vatm",sector);
                tempx_vec = reshape(tempx(sector_mask),numel(tempx(sector_mask)),1);
                tempy_vec = reshape(tempy(sector_mask),numel(tempy(sector_mask)),1);
                temp_vec = sqrt(tempx_vec.^2 + tempy_vec.^2);
                X_temp = [X_temp, temp_vec];  
            else
                temp = data_format_sector(filename,var_list{i},sector);
                temp(~ice_mask) = NaN;
                temp_vec = reshape(temp(sector_mask),numel(temp(sector_mask)),1);
                X_temp = [X_temp, temp_vec];
            end
        end 
        X_temp = [X_temp, lat_vec, lon_vec];
        [idx(j),~] = size(X_temp);
        X_total(:,:,j) = X_temp; %[X_total; X_temp];
    end
   
end


function var_out = variable_dict(var_in)
    for i = 1:numel(var_in)
        if var_in{i} == "aice"
            var_out{i} = 'SIC';
        elseif var_in{i} == "hi"
            var_out{i} = 'Ice thickness';
        elseif var_in{i} == "iage"
            var_out{i} = 'Ice age';
        elseif var_in{i} == "fsdrad"
            var_out{i} = 'Floe mean radius';
        elseif var_in{i} == "pancake"
            var_out{i} = 'Pancake ice';
        elseif var_in{i} == "large"
            var_out{i} = 'Sheet ice';
        elseif var_in{i} == "wave_sig_ht"
            var_out{i} = 'Sig. wave height';
        elseif var_in{i} == "frazil"
            var_out{i} = 'Frazil';
        elseif var_in{i} == "frzmlt"
            var_out{i} = 'Freeze/melt pot.';
        elseif var_in{i} == "congel"
            var_out{i} = 'Congelation rate';
        elseif var_in{i} == "ice_present"
            var_out{i} = 'Ice present';
        elseif var_in{i} == "peak_period"
            var_out{i} = 'Peak period';
        elseif var_in{i} == "vel"
            var_out{i} = 'Ice speed';
        elseif var_in{i} == "strint"
            var_out{i} = 'Internal stress';
        elseif var_in{i} == "Tair"
            var_out{i} = 'Air temp.';
        elseif var_in{i} == "dafsd_latg"
            var_out{i} = 'Lat growth';
        elseif var_in{i} == "dafsd_latm"
            var_out{i} = 'Lat melt';
        elseif var_in{i} == "dafsd_newi"
            var_out{i} = 'New ice';
        elseif var_in{i} == "dafsd_weld"
            var_out{i} = 'Welding';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "dafsd_wave"
            var_out{i} = 'Wave breakup';
        elseif var_in{i} == "full_fsd"
            for j = 1:12
                var_out{i+j-1} = sprintf('FSD cat. %g',j);
            end
        elseif var_in{i} == "full_fstd"
            count = 1;
           for nc = 1:5
                for j = 1:12
                    var_out{i+count-1} = sprintf('FSTD cat. (%g %g)',j,nc);
                    count = count + 1;
                end
            end
        else
            var_out{i} = var_in{i};
        end
    end
end
function [X_out,row_idx] = clearNaN(X_temp)
    [~,~,dep] = size(X_temp);
    for i = 1:dep
        % Step 0: removing NaN
        idx1 = isnan(X_temp(:,:,i));
        idx = logical(idx1);
        [len, wid] = size(X_temp(:,:,i));
        
        wid = wid - 2;
        count = 1;
        for j = 1:len
            if sum(idx(j,:)) == 0
                X_new(count,:,i) = X_temp(j,:,i);
                count  = count + 1;
            end
        end
        X_new_unstandard = X_new;
        %clear X_new_unstandard
        [row_idx(i),~] = size(X_new(:,:,i));
    end

    X_temp = [];
    for i = 1:dep
        X_temp = [X_temp; X_new(:,:,i)];
    end
    X_out = X_temp;

end




function [MIZ_width, class,MIZ] = calculate_miz_width(filename,sector,k_means,miz_def)
    aice = data_format_sector(filename,'aice',sector);
    [lat,lon,~,ulat,ulon] = grid_read('om2');
    ice_mask =  aice > 0.15;
    k_means(~ice_mask) = NaN;
    idx_miz = isnan(k_means(1,:));
    [len,~] = size(aice);
    for i = 1:len
        temp = k_means(i,~idx_miz);
        edge_class(i) = temp(end); % Take the last cell that > 0.15 SIC
    end
    idx = isnan(edge_class);
    class = mode(edge_class(~idx));
    MIZ = k_means == miz_def;%class;  
    clear MIZ_width
    
    for i = 1:len
        ulat_vec_temp = ulat(i,:);
        ulon_vec_temp = ulon(i,:);
        MIZ_vec_temp = MIZ(i,:);
        if sum(MIZ_vec_temp) == 0
            MIZ_width(i,1:2) = 0;
        else
            %south_point = find(MIZ_vec_temp, 1, 'first');
            %north_point = find(MIZ_vec_temp, 1, 'last');
            f = find(diff([0,MIZ_vec_temp,0]==1));
            p = f(1:2:end-1);  % Start indices
            y = f(2:2:end)-p;  % Consecutive ones counts
            max_y = find(y == max(y));
            p = p(max_y);
            y = y(max_y);
            south_point = p;
            north_point = p+y;
            temp = pathdist([ulat_vec_temp(south_point-1),ulat_vec_temp(north_point)],[ulon_vec_temp(south_point),ulon_vec_temp(north_point)],'km');
            MIZ_width(i,1:2) = temp(1:2);
        end
    end
    MIZ_width = MIZ_width(:,2);
end









