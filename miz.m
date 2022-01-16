%% Determining the MIZ
% The MIZ is calculated from a combination of FSD radius, significant wave
% height, and the change in FSD due to waves.
clear all
close all
addpath functions
addpath packages/bedmap2_toolbox_v4
color_map = seaicecolormap();

fsd_max = 10;
movie_miz = 0; % 0 is off, 1 is on
average_miz_switch = 0; 
sum_miz_switch = 1;
swh_ice_miz_switch = 1;
swh_ice_miz_switch2 = 0;
iceedge = 0;
wave_height = 0.01;
SIC = 0.15;
%% Preamble
user = 'noahday'; %a1724548, noahday, Noah
case_name = '12monthswim';
grid = 'gx1'; 
variable = 'fsdrad'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir, hi, uvel, vvel
time_period = 'd'; %'1','d','m','y'
season = "winter";
% Winter has 92 days
% Autumn has 92 days
% Summer has 59 (no Decemeber)
day = 1;
if season == "summer"
    datapoints = 59;
    month = 1;
elseif season == "autumn"
    datapoints = 92;
    month = 3;
elseif season == "winter"
    datapoints = 92;
    month = 6;
elseif season == "spring"
    datapoints = 91;
    month = 9;
end
year = 2005;
date = sprintf('%d-0%d-0%d', year, month, day);
map_type = 'eqaazim'; 
% Load the grid, ice shelf, and MIZ statistics.
% Grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');
dim = 2;
lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);
lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

% Load ice shelf data
shelf = bedmap2_data('icemask');
shelf = 0.0*(shelf == 1); 
shelf(shelf==0) = NaN; % 1 iceshelf
[latshelf,lonshelf] = bedmap2_data('latlon');


%% 1. CICE data in mean
if average_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month, day);
    [len, wid] = size(lat);
    % a. FSD radius < fsd_max
    dim = 2;
    variable = "fsdrad";
    fsd_data = aggregate_data(case_name,date,datapoints,variable,dim);
    idx = fsd_data > fsd_max; 
    fsd_miz = fsd_data;
    fsd_miz(idx) = 0.0;
    idx = fsd_miz > eps;
    fsd_miz(idx) = 100.0;

    % b. SWH > 0

    % dim = 2;
    % variable = "wave_sig_ht";
    % swh_data = aggregate_data(case_name,date,datapoints,variable,dim);
    % idx = swh_data > 0;
    % swh_miz = swh_data;
    % swh_miz(~idx) = 0.0;
    % idx = swh_miz > eps;
    % swh_miz(idx) = 50.0;

    % c. dafsd_wave > 0

    % dim = 3;
    % variable = "dafsd_wave";
    % wave_data = aggregate_data(case_name,date,datapoints,variable,dim);
    % idx = wave_data > 0;
    % dafsd_miz = wave_data;
    % dafsd_miz(~idx) = 0.0;
    % idx = dafsd_miz > eps;
    % dafsd_miz(idx) = 10.0;

    % Combining the datasets 
    total_miz = fsd_miz;% + swh_miz; %+ dafsd_miz;
    colour_bar = 1;
    plot_map(lat,lon,total_miz,latshelf,lonshelf,shelf,text,1,colour_bar)
end

% CICE data at each day during the season - movie
if movie_miz == 1
    figure(3)
    video_name = strcat(variable, '_', case_name, '_', '2021_11_30', '.avi');
    writerObj = VideoWriter(video_name);
    time_period = 'd'; %'1','d','m','y'
    set(writerObj,'FrameRate',datapoints/10); % 0.5 = 2 seconds per frame
    open(writerObj);

    case_name = strcat('cases/',case_name);
    plot_title_vec = " ";
    for i = 1:datapoints
       % Get the file name
       filename = strcat(case_name,"/history/iceh.",date,".nc");
       % Plot the map
       map_creator_miz(filename, plot_title_vec, i, variable, grid, map_type, user, case_name)
       date = update_date(date);
    end    

    % iterate over each image
     for k=1:datapoints
         % use imread to read the image
         figname = sprintf('frames/image%d.png',k);
         img = imread(figname);
         % convert the image to a frame using im2frame
         frame = im2frame(img);
         % write the frame to the video
         writeVideo(writerObj,frame);
     end
     % close the writer
     close(writerObj);
end % if movie

%% 2. Proportion of winter that cells satisy the FSD requirements over winter
if sum_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month, day);
    for i = 1:datapoints
       % Get the file name
       filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
       data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
       data = rearrange_matrix(data,37,dim);
       data_formatted = [data; data(end,:)];
       idx = data_formatted > fsd_max;    
       fsd_miz = data_formatted;
       fsd_miz(idx) = 0.0;
       idx = fsd_miz > eps;
       fsd_miz(idx) = 1.0;
       % only print FSD where there is significant ice
       mask = ice_mask(case_name,date,grid,SIC);
       total_data(:,:,i) = fsd_miz.*mask;
       date = update_date(date);
    end    
    sum_miz = sum(total_data,3)/datapoints;
    fsdradius = sprintf("with floe radius less than %d m", fsd_max);
    text = strcat("Proportion of days with ",fsdradius, " during ",season, " in 2005");
    colour_bar = 1;
    plot_map(lat,lon,sum_miz,latshelf,lonshelf,shelf,text,2,colour_bar)
end

%% 3. Proportion of days with SWH > wave_height within sea ice regions (greater than 15% SIC)
if swh_ice_miz_switch == 1
   date = sprintf('%d-0%d-0%d', year, month, day);
% a) Find the ice edge
    

% b) SWH > 0
    dim = 2;
    variable = "wave_sig_ht";

     for i = 1:datapoints
       % Get the file name
       filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
       data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
       data = rearrange_matrix(data,37,dim);
       data_formatted = [data; data(end,:)];
       idx = data_formatted > fsd_max;    
       swh_miz = data_formatted;
       swh_miz(idx) = 0.0;
       idx = swh_miz > wave_height;
       swh_miz(idx) = 1.0;
       % only print SWH in areas where there is sea ice
       mask = ice_mask(case_name,date,grid,SIC);
       total_data(:,:,i) = swh_miz.*mask;
       date = update_date(date);
     end
    % Proportion of days per season (datapoints)
    swh_miz = sum(total_data,3)/datapoints;
    text = "Proportion of days with SWH $ > 1$ m and SIC $> 0.15$ over winter 2005";
    colour_bar = 1;
    plot_map(lat,lon,swh_miz,latshelf,lonshelf,shelf,text,3,colour_bar)
end

%% 3b. Average SWH within sea ice regions (greater than 15% SIC)
if swh_ice_miz_switch2 == 1
   date = sprintf('%d-0%d-0%d', year, month, day);
% a) Find the ice edge
    %mask = ice_mask(case_name,date,grid,SIC);

% b) SWH > 0
    dim = 2;
    variable = "wave_sig_ht";
    swh_data = aggregate_data(case_name,date,datapoints,variable,dim);
     
    idx = swh_data < wave_height;
    swh_miz = swh_data;
    swh_miz(idx) = 0.0;
    
    % only print SWH in areas where there is sea ice
    %swh_miz = swh_miz;%.*mask;
    text = "Proportion of days with SWH $ > 1$ m and SIC $> 0.15$ over winter 2005";
    colour_bar = 1;
    plot_map(lat,lon,swh_miz,latshelf,lonshelf,shelf,text,7,colour_bar)
end

%% Ice edge
if iceedge == 1
    edge = ice_edge(case_name,date,grid);
    edge_data = zeros(size(lon));
    for i = 1:len
        edge_data(i,edge(i)) = 1.0;
    end
    figure(3)
    plot_map(lat,lon,edge_data,latshelf,lonshelf,shelf,text,3,colour_bar)
end

%% Correlation between FSD and SWH MIZ definitions
if sum_miz_switch == 1 && swh_ice_miz_switch == 1
    
    date = sprintf('%d-0%d-0%d', year, month, day);
    filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
    land_mask = ncread(filename, "tmask");
    land_mask = rearrange_matrix(land_mask,37,dim);
    combined_data = swh_miz;%~land_mask*NaN;
    correlation_data = swh_miz;
    idx = swh_miz > eps;
    combined_data(idx) = 0.0;
    correlation_data(idx) = 0.0;
    count1= 0;
    count2 = 0;

    if swh_ice_miz_switch == 1
        if sum_miz_switch == 1
            [len, wid] = size(swh_miz);
            for i = 1:len
                for j = 1:wid
                    if isnan(sum_miz(i,j)) == 0
                        if sum_miz(i,j) > eps && swh_miz(i,j) > eps
                            combined_data(i,j) = 1.0;
                        else
                            combined_data(i,j) = 0.0;
                        end
                        if sum_miz(i,j) < eps && swh_miz(i,j) < eps
                            correlation_data(i,j) = 1.0;
                            count1 = count1 + 1;
                        elseif sum_miz(i,j) > eps && swh_miz(i,j) > eps
                            correlation_data(i,j) = 1.0;
                            count2 = count2 + 1;
                        else
                            correlation_data(i,j) = 0.0;
                        end
                    end
                end
            end
            figure(4)
            text = "Combined MIZ";
            colour_bar = 1;
            plot_map(lat,lon,combined_data,latshelf,lonshelf,shelf,text,4,colour_bar)

            text = "Correlation between the two MIZ definitions";
            %plot_map(lat,lon,correlation_data,latshelf,lonshelf,shelf,text,5)
        end
    end



    % Correlation
    nan_idx = isnan(correlation_data);
    size_fsd = sum(sum(sum_miz(nan_idx)));
    size_swh = sum(sum(swh_miz(nan_idx)));

    indic = sum_miz > eps;
    indic_fsd = sum_miz;
    indic_fsd(indic) = 1.0;

    indic = swh_miz > eps;
    indic_swh = swh_miz;
    indic_swh(indic) = 2.0;

    if size_fsd > size_swh % FSD MIZ is bigger
        diff_data = sum_miz - swh_miz;
    else
        diff_data = swh_miz - sum_miz;
    end

    combined_miz = sum_miz + swh_miz;
    nan_idx = isnan(combined_miz);
    idx = combined_miz>eps;
    combined_miz(idx) = 1.0;
    total = sum(sum(combined_miz(~nan_idx)));

    add_miz = indic_fsd + indic_swh;
    idx = add_miz == 3.0;
    same_count = sum(sum(idx));
    idx = add_miz == 1.0;
    fsd_not_swh = sum(sum(idx));
    idx = add_miz == 2.0;
    swh_not_fsd = sum(sum(idx));
    %same_count = sum(sum(correlation_data(~nan_idx)));
    %total = numel(correlation_data(~nan_idx));

    correlation = same_count/total;
    fsd_prop = fsd_not_swh/total;
    swh_prop = swh_not_fsd/total;
    fprintf("The correlation between the two MIZ definitions is: %f \n",correlation);
    fprintf("Proportion that is in FSD definition and not in SWH is: %f \n",fsd_prop);
    fprintf("Proportion that is in SWH definition and not in FSD is: %f \n",swh_prop);



    % FSD - SWH

    colour_bar = 0;
    text = "Comparison of FSD and SWH MIZs; yellow: shared, green: FSD, light blue: SWH";
    plot_map(lat,lon,add_miz,latshelf,lonshelf,shelf,text,6,colour_bar)
end
%% MIZ width
 
%% Printing ice area
%   figure(4)
%  dim = 2;
%  variable = "aice";
%  aice_data = aggregate_data(case_name,date,datapoints,variable,dim);
% %addpath functions
% color_map = seaicecolormap();
% latitude = [-90,-30];
% longitude = [-180,180];
% w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [-90 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'on')
%     setm(w, 'parallellabel', 'off')
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(lat,lon,combined_data)
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     pcolorm(latshelf,lonshelf,shelf)  
%     colorbar()
%     text = strcat("The ", season, " marginal ice zone in 2005");
%     title(text,'interpreter','latex','FontSize', 18)
    
%% Functions

function ave_data = aggregate_data(case_name,date,datapoints,variable,dim)
    if dim == 2
        for i = 1:datapoints % number of days
           filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
           data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
           data = rearrange_matrix(data,37,dim);
           data_formatted = [data; data(end,:)];
           grid = "gx1";
           mask = ice_mask(case_name,date,grid,SIC);
           total_data(:,:,i) = data_formatted.*mask;
           % update date
           date = update_date(date);
        end
        ave_data = mean(total_data,3);
    else % dim == 3
        for i = 1:datapoints % number of days
            filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
            data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
            nfsd = ncread(filename, "NFSD");
            data1D = sum(data,3)/numel(nfsd); % diving by the number of categories
            data1D = rearrange_matrix(data1D,37,2);
            total_data(:,:,i) = [data1D; data1D(end,:)];
            % update date
            date = update_date(date);
        end
        ave_data = mean(total_data,3);
    end
end

function edge = ice_edge(case_name,date,grid)
    %% Find the ice edge
    % Find the ice edge and land edge
    filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
    latit = 100;
    variable = "aice";
    dim = 2;
    datapoints = 92;
    [lat,lon,row] = grid_read(grid);
    aice_data(:,:) = aggregate_data(case_name,date,datapoints,variable,dim);
    [len, wid] = size(aice_data);
    for i = 1:len
        long_ice_edge = aice_data(i,1:latit);
        pos = find(long_ice_edge < eps);
        ice_pos(i) = pos(1);
    end
    edge = ice_pos;
end


function mask = ice_mask(case_name,date,grid,conc)
    %% Find the ice edge
    % Find the ice edge and land edge
    filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
    latit = 100;
    variable = "aice";
    dim = 2;
    datapoints = 92;
    [lat,lon,row] = grid_read(grid);
    filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
    aice_data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
    aice_data = rearrange_matrix(aice_data,37,dim);
    aice_data_formatted = [aice_data; aice_data(end,:)];
    [len, wid] = size(aice_data_formatted);
    ice_pos = zeros(len,wid);
    for i = 1:len
        long_ice_edge = aice_data_formatted(i,1:latit);
        pos = find(long_ice_edge > conc);
        ice_pos(i,pos) = 1;
    end
    mask = ice_pos;
end

function y = map_creator_miz(filename, plot_title_vec, i, variable, grid, map_type, user, case_name)
%MAP_CREATOR creates map images of Antarctica given netcdf data files
%   filename: the directory to the .nc file
%   plot_title: string containing the title for the plot
%   i: the date of the data file
%   variable: the variable in the dataset we want to plot
%   grid: specify the grid. eg. 'gx1', 'gx3'

% Load ice shelf data
addpath packages/bedmap2_toolbox_v4

shelf = bedmap2_data('icemask');
shelf = 0.0*(shelf == 1); 
shelf(shelf==0) = NaN; % 1 iceshelf
[latshelf,lonshelf] = bedmap2_data('latlon');

% grid
    dim = 2;
if grid == 'gx3'
    row = 11;

    ulat = ncread('grid_gx3.nc','ulat');
    ulon = ncread('grid_gx3.nc','ulon');

    % converting to degrees
    lon = rad2deg(ulon);
    lat = rad2deg(ulat);

    lat = rearrange_matrix(lat,row,dim);
    lon = rearrange_matrix(lon,row,dim);
   
else
    row = 37;
    lat = ncread('grid/global_gx1.bathy.nc','TLAT');
    lon = ncread('grid/global_gx1.bathy.nc','TLON');

    lat = rearrange_matrix(lat,row,dim);
    lon = rearrange_matrix(lon,row,dim);


    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
        
end

%filename = strcat('cases/',filename);
data = ncread(filename, variable);
[~, ~, n] = size(data);

for level = 1:n
    data_1 = data(:,:,level);

    latitude = [-90,90];
    longitude = [-180,180];

    data_1 = rearrange_matrix(data_1,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_1 = [data_1; data_1(end,:)];
    %% Threshold for the data
    variable = "fsdrad";
    fsd_max = 50.0;
    idx = data_1 > fsd_max; 
    fsd_miz = data_1;
    fsd_miz(idx) = 0.0;
    idx = fsd_miz > eps;
    fsd_miz(idx) = fsd_max;    
    %% Mapping
    color_map = seaicecolormap();
if map_type == 'cassini'
    x_origin = -20;
else
    x_origin = -90;
end
    w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [x_origin 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,fsd_miz)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    pcolorm(latshelf,lonshelf,shelf)  

 if variable == "fsdrad"
         plot_variable = "FSD radius ";
         unit = "metres";
     elseif  variable == "fsdrad_d"
         plot_variable = "FSD radius ";
         unit = "metres";
     elseif variable == "wave_sig_ht"
         plot_variable = "Significant wave height ";
         unit = "metres";
     elseif variable == "wave_sig_ht_d"
         plot_variable = "Significant wave height ";
         unit = "metres";
     elseif variable == "peak_period"
         plot_variable = "Peak period ";
         unit = "s";
     elseif variable == "peak_period_d"
         plot_variable = "Peak period ";
         unit = "s";
     elseif variable == "aice"
         plot_variable = "Concentration of ice ";
         unit = " ";
     elseif variable == "aice_d"
         plot_variable = "Concentration of ice ";
         unit = " ";
     elseif variable == "mean_wave_dir"
         plot_variable = "Mean wave direction (rads) ";
         unit = "radians";
     elseif variable == "mean_wave_dir_d"
         plot_variable = "Mean wave direction (rads) ";
         unit = "radians";
     else
         plot_variable = variable;
 end
    set(gcf, 'Position',  [100, 100, 1000, 800])
    set(gcf,'Visible', 'off')
    fontSize = 20; 
    plot_title = strcat(plot_variable, plot_title_vec);
    title(plot_title, 'FontSize', fontSize);
    a=colorbar;
    label_c = ylabel(a,unit,'FontSize',16,'Rotation',270);
    label_c.Position(1) = 4;
    label_h.Position(2) = 1; % change vertical position of ylabel
    limit = colorlims(variable);
    caxis(limit);
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    saveas(gcf,fullfile(filedir, figname));
end
end

function [] = plot_map(lat,lon,total_miz,latshelf,lonshelf,shelf,text,i,colourbar)
    %% Mapping
    latitude = [-90,-30];
    longitude = [-180,180];
    addpath functions
    figure(i)
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-30]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 30);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(latshelf,lonshelf,shelf)
        pcolorm(lat,lon,total_miz)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])  
        colorbar
        if colourbar == 1
            caxis([0,1])
        end
        title(text,'interpreter','latex','FontSize', 18)
        set(gcf, 'Position',  [100, 100, 1000, 800])
end