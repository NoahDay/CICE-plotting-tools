% Classification of CICE ice covers

% We want to experiment with different combination of variables over
% multiple years of data. 
% 1. Only ice cover "snapshot" variables
% 2. Add dynamic terms
% 3. Add thermodynamic terms
% 4. All of the terms

% Each case will have its own method:
% 1. Get filenames
% 2. Read in data
% 3. Clear NaNs
% 4. Standardise the data
% 5.Classify using kmeans

% Setup
clear all
clc
close all

%%
addpath functions
addpath /Users/noahday/GitHub/CICE-analyser/processing

historydir = '/Volumes/NoahDay5TB/WIMonAlessandroRun/history/2019/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end


%% Select the relevant variables
% Columns are data
sector = "SH";
close all
var_list = {'aice','hi','hs','sice','ice_present','iage','alvl','vlvl','ardg','vrdg','FYarea'};
[X_raw, row_idx]= read_data_vec(filenames,sector,var_list); % [var_list, lon, lat]
%[~, ~,dep] = size(X_temp);

%clear X_temp
label_vec = variable_dict(var_list);
size(X_raw)
%
data.Xunstandard = X_raw;
data.row_idx = row_idx;
save_filename = strcat('coversnapshot_5percent_2019.mat');
save(save_filename,'data','-v7.3');
clear data

%% Load in data
clear X row_idx
load('X_waves.mat')
X.waves = data.Xunstandard;
row_idx.waves = data.row_idx;
clear data


%% Clean the data
% NEED TO DO THIS ALL TOGETHER SO ITS ALL STANDARDISED CORRECTLY

[Xnan,row_idx] = clearNaN(X_raw);
disp('Clear NaN done!')

dimension = size(Xnan);
% Xnan = NaN's have been removed.
%%
%Xnan.all = Xnan.waves;
%%
cleaned_data.Xnan = Xnan;
cleaned_data.label_vec = label_vec;
cleaned_data.row_idx = row_idx;
cleaned_data.dimension = dimension;
save('cleanded_cover_data2019.mat','cleaned_data','-v7.3');

%%
load('cleanded_kmeans_data.mat')
Xnan.all = cleaned_data.X_no_nan;
dimension = cleaned_data.dimension;
row_idx = cleaned_data.row_idx;
clear cleaned_data
%% Standardise the data
X_standard_all = Xnan;
[~,wid] = size(X_standard_all);

for j = 1:wid-2 % Don't standardise latitude and longitude
    if min(X_standard_all(:,j)) < 0
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
disp('Standardisation done!')
%clear X_temp
%%
eva = evalclusters(X_standard_all(:,1:end-2),'kmeans','CalinskiHarabasz','KList',1:5);
temp = eva.CriterionValues;
conFigure(11)
f = figure;
plot(eva)
exportgraphics(f,'calinskiHarabasz_cover.pdf','ContentType','vector')

%% k-means clustering
rng(2022)
X = X_standard_all(:,1:end-2);
num_clusters = 3;
tic
[idx,C] = kmeans(X,num_clusters,'MaxIter',300);
toc
% %%
% kmeans_cluster.idx = idx;
% kmeans_cluster.row_idx = row_idx;
% kmeans_cluster.label_vec = label_vec;
% kmeans_cluster.C = C;
% kmeans_cluster.X_standard_all = X_standard_all;
% kmeans_cluster.num_clusters = num_clusters;
% 
% save_filename = strcat('kmeans_3cases_',sprintf('%g',num_clusters),'_classes.mat');
% save(save_filename,'kmeans_cluster','-v7.3');
% %%
% %save_filename = strcat('kmeans_3cases_',sprintf('%g',num_clusters),'_classes.mat');
% %load(save_filename);
% 
% C = kmeans_cluster.C;
% idx = kmeans_cluster.idx;
% row_idx = kmeans_cluster.row_idx;
% X_standard_all = kmeans_cluster.X_new_unstandard;%X_standard_all;
% num_clusters = kmeans_cluster.num_clusters;
% label_vec = kmeans_cluster.label_vec;
%% Average stats
average_stats = [];
X_temp = [];
[~,~,dep] = size(Xnan);
X_temp = Xnan;
X_temp = [X_temp(:,1:end-2), idx];

for i = 1:num_clusters
    index_class_data = X_temp(:,end) == i;
    class_data = X_temp(index_class_data,1:end-1);
    average_stats(i,:) = mean(class_data);
end

close all
conFigure(12,10)
f = figure('Position',[0, 0, 60, 10]);
for i = 1:length(label_vec)
  subplot(2,ceil(length(label_vec)/2),i)
  bar(average_stats(:,i))
  xticks(1:length(average_stats(:,i)))
  title(label_vec{i})
end
exportgraphics(f,strcat('stat_comparison_',sprintf('%g',num_clusters),'_cover_5_percent_clusters.pdf'),'ContentType','vector')

%% FSD plots
close all
NFSD = ncread(filenames.waves(1,:),"NFSD");
NCAT = ncread(filenames.waves(1,:),"NCAT");
N = num_clusters;
for i = 1:num_clusters
    fstd(1,:,i) = average_stats(i,15:26);
    fstd(2,:,i) = average_stats(i,27:27+11);
    fstd(3,:,i) = average_stats(i,39:39+11);
    fstd(4,:,i) = average_stats(i,51:51+11);
    fstd(5,:,i) = average_stats(i,63:63+11);
end
clear fsd
for i = 1:num_clusters
    fsd(i,:) = sum(fstd(:,:,i));
end

region_label = [{'1'},{'2'},{'3'},{'4'},{'5'}];
conFigure(11)
f = figure;
N = 5;
C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
Ctemp(5,:) = C1(5,:);
Cmap = Ctemp([1,2,4,3,5],:);
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
set(gcf,'Visible', 'on')
loglog(NFSD,fsd([1],:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label,'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
ylim([10^-5,10^2])
%exportgraphics(f,'floesizecomp1.pdf','ContentType','vector')

axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
set(gcf,'Visible', 'on')
loglog(NFSD,fsd([1,2],:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label,'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
ylim([10^-5,10^2])
%exportgraphics(f,'floesizecomp2.pdf','ContentType','vector')

axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
set(gcf,'Visible', 'on')
loglog(NFSD,fsd([1,2,3],:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label,'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
ylim([10^-5,10^2])
%exportgraphics(f,'floesizecomp3.pdf','ContentType','vector')

axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
set(gcf,'Visible', 'on')
loglog(NFSD,fsd(:,:),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label,'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
ylim([10^-5,10^2])
%exportgraphics(f,'floesizecomp4.pdf','ContentType','vector')

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

index_lat = X_new(:,end-1) == X_new(30,end-1);
index_lon = X_new(:,end) == X_new(30,end);
index_both = index_lat.*index_lon;
sum(index_both)
%% Worldmap
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;
num_clusters = 5;
%dimension.pancake(1) = cleaned_data.dimension.pancake(1);
%dimension.sheet(1) = cleaned_data.dimension.sheet(1);

addpath functions
[lat,lon,~,ulat,ulon] = grid_read('om2');
clear X_map MIZ_width
% Make worldmap with colour matching
uarea = data_format(filenames.waves(1,:),"uarea");

[num_files,~] = size(filenames.waves);
sector = "SH";
coords = sector_coords(sector);
font_size = 5;
plot_type = "kmeans";
plotting = "on";
conFigure(11)
%X_standard_all(:,end-1:end) = Xnan.all(:,end-1:end);
%X_new = X_standard_all(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1),:);
%X_new = Xnan.waves;
X_new = Xnan.all(1:dimension.waves(1),:);

[len,wid] = size(lat);
N = 5;

C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
%Cmap = Ctemp([1,2,4,3],:);

file_number = 250;


index_lat = X_new(:,end-1) == X_new(1,end-1);
index_lon = X_new(:,end) == X_new(1,end);
index_both = index_lat.*index_lon;
diff_index = diff(find(index_both)');
row_file = find(index_both);%[0,cumsum(row_idx)];
row_vec = row_file(file_number)+1:row_file(file_number+1);
%temp_idx = idx(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1)); % waves:waves+sheet
temp_idx = idx(1:dimension.waves(1));
file_idx = temp_idx(row_vec);
[aice, sector_mask] = data_format_sector(filenames.waves(file_number,:),'aice',sector);
%[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
k_means = NaN.*ones(len,wid);
for i = 1:length(file_idx)
    [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
    k_means(lon_pos,lat_pos) = file_idx(i); 
    k_means(~sector_mask) = NaN;
end

ice_mask =  aice > 0.15;
% k_means(~ice_mask) = NaN;
if plotting == "on"
    f = figure;
    set(gcf,'Visible', 'on')
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'grid', 'off');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,k_means)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        %colorbar; %cmocean('deep');
        colormap(Cmap)
        cb = colorbar;
        title(dirdates.waves(file_number,:),'Color','black','FontSize',font_size+5)
        %plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
        
        if plot_type == "pc1"
            caxis([0,1])
        elseif plot_type == "kmeans"
            %cb = colorbar; cmocean('deep',num_clusters)
            %cb.TickLabels = region_label; 
            cb.Ticks = 1:num_clusters;
            %cb.Location = 'southoutside';%cb.Location = 'northoutside';
            cb.Location = 'eastoutside';
            caxis([1,num_clusters])
            cb.AxisLocation = 'in';
        end
end


%% MOVIES
% [waves; sheet; pancake]
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;
num_clusters = 4;
%dimension.pancake(1) = cleaned_data.dimension.pancake(1);
%dimension.sheet(1) = cleaned_data.dimension.sheet(1);

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
plotting = "off";
conFigure(11)
%X_standard_all(:,end-1:end) = Xnan.all(:,end-1:end);
%X_new = X_standard_all(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1),:);
X_new = Xnan.all(1:dimension.waves(1),:);
%X_new = Xnan.all(1:dimension.pancake(1),:);

[len,wid] = size(lat);
N = num_clusters;

C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
Cmap = Ctemp([1,2,4,3],:);


file_number_vec = 1:1:365;
for file_number = file_number_vec %num_files-1
    %row_file = [0,cumsum(row_idx.waves)];
    if plot_type == "pc1"
        file_idx = score(row_file(file_number)+1:row_file(file_number+1),1);
        file_idx = file_idx + abs(min(file_idx));
        file_idx = file_idx./max(file_idx);
    elseif plot_type == "kmeans"
        index_lat = X_new(:,end-1) == X_new(30,end-1);
        index_lon = X_new(:,end) == X_new(30,end);
        index_both = index_lat.*index_lon;
        diff_index = diff(find(index_both)');
        %row_file = [0,cumsum(diff_index)];
        row_file = find(index_both);%[0,cumsum(row_idx)];
        row_vec = row_file(file_number)+1:row_file(file_number+1);
        %temp_idx = idx(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1)); % waves:waves+sheet
        temp_idx = idx(1:dimension.waves(1));
        file_idx = temp_idx(row_vec);
    end
    [aice, sector_mask] = data_format_sector(filenames.waves(file_number,:),'aice',sector);
    %[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    
    ice_mask =  aice > 0.15;
   % k_means(~ice_mask) = NaN;
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
            setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
            setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
            setm(w, 'mlabelparallel', 'south','FontColor','black','FontSize',3);
            setm(w, 'grid', 'off');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,k_means)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
            %colorbar; %cmocean('deep');
            colormap(Cmap)
            cb = colorbar;
            title(dirdates.waves(file_number,:),'Color','black','FontSize',font_size+5)
            %plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
            
            if plot_type == "pc1"
                caxis([0,1])
            elseif plot_type == "kmeans"
                %cb = colorbar; cmocean('deep',num_clusters)
                %cb.TickLabels = region_label; 
                cb.Ticks = 1:num_clusters;
                %cb.Location = 'southoutside';%cb.Location = 'northoutside';
                cb.Location = 'eastoutside';
                caxis([1,num_clusters])
                cb.AxisLocation = 'in';
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

%% Save the files
ap_pancake.perimeter_region = perimeter_region;
ap_pancake.area_region = area_region;
ap_pancake.xvec = file_number_vec;
ap_pancake.MIZ_width = MIZ_width;
save_filename = strcat('a_p_pancake.mat');
save(save_filename,'ap_pancake','-v7.3');



%%
% % cbar = cmocean('deep',num_clusters);
% % close all
% % N = 6;
% % %cmocean()
% % file_number_vec = 1:351;
% % region_label = {'1','2','4'};
% % date_label = datetime(dirdates.waves);
% % date_label_vec = date_label(file_number_vec);
% % 
% % 
% % conFigure(11,0.2)
% % f = figure;
% % subplot(3,1,2)
% % C = linspecer(N);
% % axes('NextPlot','replacechildren', 'ColorOrder',C);
% % plot(date_label_vec,[perimeter_region([1,2,4],file_number_vec); sum(perimeter_region(:,file_number_vec))],'LineWidth',5);
% % %set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
% % ylabel('Perimeter [block edges]')
% % %legend(region_label,'Location','southeast')
% % 
% % subplot(3,1,3)
% % plot(date_label_vec,[area_region([1,2,4],file_number_vec).*1e-6.*1e-6; sum(area_region(:,file_number_vec).*1e-6.*1e-6)],'LineWidth',5)
% % ylabel('Area [million km$^2$]')
% % %legend(region_label,'Location','southeast')
% % 
% % % Create a tile on the right column to get its position
% % legend_label = region_label;
% % legend_label{end+1} = 'Sum of regions';
% % ax = subplot(3,1,1,'Visible','off');
% % axPos = ax.Position;
% % delete(ax)
% % % Construct a Legend with the data from the sub-plots
% % hL = legend(legend_label);
% % % Move the legend to the position of the extra axes
% % hL.Position(1:2) = axPos(1:2);






% subplot(4,1,4)
% perim_per_area = perimeter_region([1,2,4],file_number_vec)./(area_region([1,2,4],file_number_vec).*1e-6.*1e-6);
% plot(date_label_vec,perim_per_area,'LineWidth',5)
% ylabel('Perimeter per Area [block per million km$^2$]')



%text(1,1,sprintf("$sigma =$ %g",std(perimeter_region'./(area_region./1000)')))
% ta = annotation('textarrow', [0.76 0.83], [0.71 0.71]);
% ta.String = sprintf("$\sigma =$ %g",std(perimeter_region(4,:)'./(area_region(4,:)./1000)'));
% ta.Interpreter = "latex";
%legend(region_label,'Location','southeast')



% diff_perim = perim_per_area(:,2:end)' - perim_per_area(:,1:end-1)';
% subplot(2,3,5)
% bar(std(diff_perim))
% xticklabels(region_label)
% ylabel('S.D. $\Delta$ perimeter per area')



%exportgraphics(f,'wave_area_perimeter.pdf','ContentType','vector')

%%
% LINE COLORS
close all
N=6;
X = linspace(0,pi*3,1000);
Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N);

subplot(2,2,1);
nexttile
C = linspecer(N);
axes('NextPlot','replacechildren', 'ColorOrder',C);
p = plot(date_label_vec,[perimeter_region([1,2,4],file_number_vec); sum(perimeter_region(:,file_number_vec))],'LineWidth',5);
%ylim([-1.1 1.1]);
%%
%file_number_vec = 1:361;
% load('a_p_sheet.mat');
% area_region = ap_sheet.area_region;
file_number_vec = 1:364;
date_label = datetime(dirdates.waves);
date_label_vec = date_label(file_number_vec);
close all
conFigure(11,1.4)
% f = figure;
N = 4;
C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
Cmap = Ctemp([1,2,4,3],:);
% axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
% plot(date_label_vec,[perimeter_region([1,2,3,4],file_number_vec); ],'LineWidth',3,'LineStyle','--');
% %set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
% ylabel('Perimeter [block edges]')
%legend(region_label,'Location','southeast')

%exportgraphics(f,'wave_perimeter.pdf','ContentType','vector')

f = figure;
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(date_label_vec,[area_region([1,2,3,4],file_number_vec).*1e-6.*1e-6],'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)
%set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')

%exportgraphics(f,'wave_area.pdf','ContentType','vector')
%%
file_number_vec = 1:364;
Cmap(1,:) = C1(2,:);
Cmap(2,:) = C1(1,:);
f = figure;
axes('NextPlot','replacechildren', 'ColorOrder',Cmap);
plot(date_label_vec,[area_region([4],file_number_vec).*1e-6.*1e-6; sum(area_region([1,2],file_number_vec)).*1e-6.*1e-6],'LineWidth',3,'LineStyle','-') % sum(area_region(:,file_number_vec).*1e-6.*1e-6)
%set(p, {'color'}, {[cbar(1,:)];[cbar(2,:)];[cbar(4,:)];[0.5,0.1,0.1]});
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')

exportgraphics(f,'pack.pdf','ContentType','vector')

%%
subplot(2,1,2)
C = linspecer(N);
axes('NextPlot','replacechildren', 'ColorOrder',C);
plot(date_label_vec,[area_region([1,2,4],file_number_vec).*1e-6.*1e-6; sum(area_region(:,file_number_vec).*1e-6.*1e-6)],'LineWidth',5)
ylabel('Area [million km$^2$]')
%legend(region_label,'Location','southeast')




% % Create a tile on the right column to get its position
% legend_label = region_label;
% legend_label{end+1} = 'Sum of regions';
% ax = subplot(3,1,1,'Visible','off');
% axPos = ax.Position;
% delete(ax)
% % Construct a Legend with the data from the sub-plots
% hL = legend(legend_label);
% % Move the legend to the position of the extra axes
% hL.Position(1:2) = axPos(1:2);




%%
close
t  = 0:5;
x1 = zeros(size(t));
x2 = zeros(size(t));
tt  = linspace(0, 5, 20);
x1t = interp1(t, x1, tt);
x2t = interp1(t, x2, tt);
plot(tt, x1t, '--o', 'MarkerIndices', 1:5:length(tt));
plot(tt, x2t, '--o', 'MarkerIndices', 3:5:length(tt));
%%


NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");
f = figure;
set(gcf,'Visible', 'on')
loglog(NFSD,average_stats([1,2,4],14:end),'LineWidth',3,'LineStyle','-','Marker','s','MarkerSize',5)
legend(region_label{[1,2,4]},'location','southwest')
ylabel('Sea ice concentration [$\%$]')
xlabel('Floe size [m]')
exportgraphics(f,'floesizecompminus3.pdf','ContentType','vector')









%% Evolution of statistics in each class
% [waves; sheet; pancake]
SIC = 0.15;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;
num_clusters = 4;
%dimension.pancake(1) = cleaned_data.dimension.pancake(1);
%dimension.sheet(1) = cleaned_data.dimension.sheet(1);

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
plotting = "off";
conFigure(11)
%X_standard_all(:,end-1:end) = Xnan.all(:,end-1:end);
%X_new = X_standard_all(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1),:);
X_new = Xnan.all(1:dimension.waves(1),:);
%X_new = Xnan.all(1:dimension.pancake(1),:);

[len,wid] = size(lat);
N = num_clusters;

C1 = linspecer(N);
Ctemp(1,:) = C1(4,:);
Ctemp(2,:) = C1(3,:);
Ctemp(3,:) = C1(1,:);
Ctemp(4,:) = C1(2,:);
Cmap = Ctemp([1,2,4,3],:);


file_number_vec = 1:30:365;
for file_number = file_number_vec
    if plot_type == "pc1"
        file_idx = score(row_file(file_number)+1:row_file(file_number+1),1);
        file_idx = file_idx + abs(min(file_idx));
        file_idx = file_idx./max(file_idx);
    elseif plot_type == "kmeans"
        index_lat = X_new(:,end-1) == X_new(30,end-1);
        index_lon = X_new(:,end) == X_new(30,end);
        index_both = index_lat.*index_lon;
        diff_index = diff(find(index_both)');
        %row_file = [0,cumsum(diff_index)];
        row_file = find(index_both);%[0,cumsum(row_idx)];
        row_vec = row_file(file_number)+1:row_file(file_number+1);
        %temp_idx = idx(dimension.waves(1)+dimension.sheet(1):dimension.waves(1)+dimension.sheet(1)+dimension.pancake(1)); % waves:waves+sheet
        temp_idx = idx(1:dimension.waves(1));
        file_idx = temp_idx(row_vec);
    end
    [aice, sector_mask] = data_format_sector(filenames.waves(file_number,:),'aice',sector);
    %[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    
    ice_mask =  aice > 0.15;
   % k_means(~ice_mask) = NaN;

    % Calculate the width of the MIZ
    % CHECK WHAT number MIZ is !!!!!!!!
    %[MIZ_width(file_number,:), miz_class, MIZ_zone] = calculate_miz_width(convertStringsToChars(filenames.waves(file_number,:)),sector,k_means,1);
    

    % Calculate evolution of each statistic
    for rgn_number = 1:num_clusters
        region_data = k_means == rgn_number;
        [len, wid] = size(k_means); % Size of grid
        wave_sig_ht_temp = data_format_sector(filenames.waves(file_number,:),'wave_sig_ht',sector);
        wave_sig_ht(rgn_number, file_number) = mean(mean(wave_sig_ht_temp.*region_data));
        ppd_temp = data_format_sector(filenames.waves(file_number,:),'peak_period',sector);
        ppd(rgn_number, file_number) = mean(mean(wave_sig_ht_temp.*region_data));
    end

    clear aice ice_mask k_means file_idx X_map
end

































%% FUNCTIONS
function [X_total, idx] = read_data_vec(filenames,sector,var_list)
    [nfiles,~] = size(filenames);
    if nfiles == 1
        filename = filenames;
    else
        filename = filenames(1,:);
    end
    
    SIC_threshold = 0.05;
    [~, sector_mask] = data_format_sector(filename,"aice",sector);
    [lat,lon,~,~,~] = grid_read('om2');
    lat_vec = reshape(lat(sector_mask),numel(lat(sector_mask)),1);
    lon_vec = reshape(lon(sector_mask),numel(lon(sector_mask)),1);

    NFSD = ncread(filename,"NFSD");
    NCAT = ncread(filename,"NCAT");
    [floe_binwidth, ~, ~, ~] = cice_parameters(NFSD);
    
    X_total = [];
    for j = 1:nfiles % Each file
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
        elseif var_in{i} == "sice"
            var_out{i} = 'Ice salinity';
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
        [len1,wid1] = size(temp);
        if len1 == 0
            edge_class(i) = 0;
        elseif wid1 == 0
            edge_class(i) = 0;
        else
            edge_class(i) = temp(end); % Take the last cell that > 0.15 SIC
        end
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









