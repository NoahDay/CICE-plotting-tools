%% Determining the MIZ
% The MIZ is calculated from a combination of FSD radius, significant wave
% height, and the change in FSD due to waves.
clear all
close all
addpath functions
addpath packages/bedmap2_toolbox_v4
% Switches
sic_miz_switch = 0;
swh_miz_switch = 0;
fsd_miz_switch = 0;
iceedge = 0;

plotting = 0;
printing = 0;
sector = [-65, 360-132.2];
ssd = 0; % ssd or local data
SIC = [0.15, 0.8];
SWH = 0.0001;
fsd_max = 50;
%% Preamble
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'twoyearproper';
if ssd == 1
    ssd_dir = '/Volumes/Noah_SSD/run_data';
    filedir = strcat(ssd_dir,case_name);
else
    filedir = strcat('cases/',case_name);
end
grid = 'gx1'; 
time_period = 'd'; %'1','d','m','y'
season = "winter";
% Winter has 92 days
% Autumn has 92 days
% Summer has 59 (no Decemeber)
day = 1;
if season == "summer"
    datapoints = 59;
    month_init = 1;
elseif season == "autumn"
    datapoints = 92;
    month_init = 3;
elseif season == "winter"
    datapoints = 92;
    month_init = 6;
elseif season == "spring"
    datapoints = 91;
    month_init = 9;
end

day = 1;
month_init = 9;
year = 2006;
date = sprintf('%d-0%d-0%d', year, month_init, day);
months = [2, 5, 9, 12];
dim = 2;
% Load the grid, ice shelf, and MIZ statistics.
% Grid
[lat,lon,row] = grid_read(grid);

% % Load ice shelf data
% shelf = bedmap2_data('icemask');
% shelf = 0.0*(shelf == 1); 
% shelf(shelf==0) = NaN; % 1 iceshelf
% [latshelf,lonshelf] = bedmap2_data('latlon');
figcount = 0;


%% 1. SIC defintion
% The SIC MIZ is defined between 0.15 and 0.8 SIC
if sic_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 2;
    variable = "aice";
    sic_data = data_format_sector(filename,variable,"SA");
    idx = sic_data < SIC(1); 
    sic_miz = sic_data;
    sic_miz(idx) = 0.0;
    idx = sic_data > SIC(2); 
    sic_miz(idx) = 0.0;
    figcount = figcount + 1;
    figure(figcount)
    map_plot(sic_data,variable,"SA",grid);
end

%% 2. SWH definition
if swh_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 2;
    variable = "wave_sig_ht";
    swh_data = data_format_sector(filename,variable,sector);
    idx = swh_data < eps; 
    swh_miz = swh_data;
    swh_miz(idx) = 0.0;
    figcount = figcount + 1;
    figure(figcount)
    map_plot(swh_miz,variable,sector,grid);
end

%% 3. FSD definition
if fsd_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 2;
    variable = "fsdrad";
    fsd_data = data_format_sector(filename,variable,sector);
    idx = fsd_data > 10; 
    fsd_miz = fsd_data;
    fsd_miz(idx) = 0.0;
    figcount = figcount + 1;
    figure(figcount)
    map_plot(fsd_miz,variable,sector,grid,[0,10]);
end
%% 4. MIZ widths
Data = struct('Month',{},'SIC',{},'SWH',{},'FSD',{});
for k = 1:4
    month_init = months(k);
    if month_init < 10
        date = sprintf('%d-0%d-0%d', year, month_init, day);
    else
        date = sprintf('%d-%d-0%d', year, month_init, day);
    end
    if month_init == 1 || month_init == 3 || month_init == 5 || month_init == 7 || month_init == 8 || month_init == 10 ||  month_init == 12
        datapoints = 31;
    elseif month_init == 2
        datapoints = 28;
    else
        datapoints = 30;
    end
    distSIC = zeros(1,datapoints);
    distSWH = zeros(1,datapoints); 
    distFSD = zeros(1,datapoints);
    for j = 1:datapoints
        filename = strcat(filedir,"/history/iceh.",date,".nc");
        % 4. a) SIC
        [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
        variable = "aice";
        data = data_format(filename,variable,row,lat,lon,dim);
        if plotting == 1
            figcount = figcount + 1;
            figure(figcount)
            [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                title('FSDrad along transect')
        else
            [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        end
        transect_data = output_data(lon_out,:); % Take the southern hemisphere
        
        
        idx1 = transect_data > SIC(1);
        idx2 = transect_data < SIC(2);
        idx_transect = logical(idx1.*idx2);
        
        points = lat(lon_out,idx_transect);
        southern_hemi_points = points(points < 0);
        [~,wid] = size(southern_hemi_points);
        if wid == 0 % No MIZ
            distSWH(j) = 0;
        else
            distSIC(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
        end
        if printing == 1
            fprintf('The width of the SIC MIZ along %g E is %g km\n', sector(2), distSIC)
        end
        
        % 4.b SWH
        [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
        variable = "wave_sig_ht";
        data = data_format(filename,variable,row,lat,lon,dim);
        if plotting == 1
            figcount = figcount + 1;
            figure(figcount)
            [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        else
            [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        end
    
        transect_data = output_data(lon_out,:); % Take the southern hemisphere
        
        idx_transect = transect_data > SWH;
        
        points = lat(lon_out,idx_transect);
        southern_hemi_points = points(points < 0);
        [~,wid] = size(southern_hemi_points);
        if wid == 0 % No MIZ
            distSWH(j) = 0;
        else
            distSWH(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
        end
        if printing == 1
            fprintf('The width of the SWH MIZ along %g E is %g km\n', sector(2), distSWH)
        end
        % 4.c FSD
        
        [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
        variable = "fsdrad";
        data = data_format(filename,variable,row,lat,lon,dim);
        if plotting == 1
            figcount = figcount + 1;
            figure(figcount)
            [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                title('FSDrad along transect')
        else
            [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
    
        end
        transect_data = output_data(lon_out,:); % Take the southern hemisphere
        idx1 = transect_data < fsd_max;
        idx2 = transect_data > eps;
        idx_transect = logical(idx1.*idx2);
        
        
        points = lat(lon_out,idx_transect);
        southern_hemi_points = points(points < 0);
        [~,wid] = size(southern_hemi_points);
        if wid == 0 % No MIZ
            distFSD(j) = 0;
        else
            distFSD(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
        end
        if printing == 1
            fprintf('The width of the FSD MIZ along %g E is %g km\n', sector(2), distFSD)
        end
        date = update_date(date);
    end
    Data(k).Month = [distSIC;distSWH;distFSD];
end
%% Plotting
%% a) Boxplot widths
f = figure;
t = tiledlayout(1,length(months));
for i = 1:length(months)
    plot_data = Data(i).Month;
    figcount = figcount + 1;
    %figure(figcount)
    nexttile
    boxplot(plot_data',["SIC","SWH","FSD"])
    ylabel('MIZ width (km)')
    xlabel('MIZ definition')
    ylim([0,1200])
    text = strcat("MIZ widths of the %g E transect over " ,monthName(months(i))," %g");
    title(sprintf(text,sector(2),year))
    f.Position = [100 100 1500 400];
end
%% b) Correlations
f = figure(2);
t2 = tiledlayout(1,length(months));
for i = 1:length(months)
    plot_data = Data(i).Month;
    figcount = figcount + 1;
    %figure(figcount)
    nexttile
    x = plot_data(2,:);
    y = plot_data(3,:);
    %scatter(x,y)
    % Get coefficients of a line fit through the data.
    coefficients = polyfit(x, y, 1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(x), max(x), 1000);
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients , xFit);
    % Plot everything.
    plot(x, y, 'b.', 'MarkerSize', 15); % Plot training data.
    hold on; % Set hold on so the next plot does not blow away the one we just drew.
    plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
    grid on;
    ylabel('SWH MIZ width (km)')
    xlabel('FSD MIZ width (km)')
    ylim([0,800])
    xlim([0,1000])
    text = strcat(monthName(months(i))," %g along %g E transect");
    title(sprintf(text,year,sector(2)))
    f.Position = [100 100 1500 400];
end
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
function name = monthName(num)
    name = month(datetime(1,num,1), 'name');
    name = name{1};
end