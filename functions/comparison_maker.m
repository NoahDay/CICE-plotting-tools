function y = comparison_maker(filename, filename2, plot_title1, plot_title2, variable, grid, map_type, user)
%MAP_CREATOR creates map images of Antarctica given netcdf data files
%   filename: the directory to the .nc file
%   plot_title: string containing the title for the plot
%   i: the date of the data file
%   variable: the variable in the dataset we want to plot
%   grid: specify the grid. eg. 'gx1', 'gx3'
dim = 2;
% grid
if grid == "gx3"
    row = 11;
    
    ulat = ncread('grid/grid_gx3.nc','ulat');
    ulon = ncread('grid/grid_gx3.nc','ulon');

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

    %% Data 1

data = ncread(filename, variable);
[~, ~, n] = size(data);

%for level = 1:n
    data_1 = data(:,:);

    latitude = [-90,90];
    longitude = [-180,180];

    data_1 = rearrange_matrix(data_1,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_1 = [data_1; data_1(end,:)];

    %% Data 2
    data2 = ncread(filename2, variable);
    [~, ~, n] = size(data2);

%for level = 1:n
    data_2 = data2(:,:);

    latitude = [-90,90];
    longitude = [-180,180];

    data_2 = rearrange_matrix(data_2,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_2 = [data_2; data_2(end,:)];
    
    %% Mapping
    
    tiledlayout(1,2) % Requires R2019b or later

% Left plot
    nexttile
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_1)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

    %antarctica = shaperead('landareas', 'UseGeoCoords', true,...
    %  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

  
    
    title(plot_title1);
    %caxis([-10 10]) %fsdrad [80 250]
    colorbar
    
% Right plot
nexttile
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_2)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

    %antarctica = shaperead('landareas', 'UseGeoCoords', true,...
    %  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

  
    
    title(plot_title2);
    %caxis([-10 10]) %fsdrad [80 250]
    colorbar
    
    
    
    %patchm(antarctica.Lat, antarctica.Lon, [1 1 1]);
    %figname = sprintf('image%d.png', i); 
    %text = sprintf('/Users/%s/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames', user);
    %fname = '/Users/noahday/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames';
    %saveas(gcf,fullfile(fname, figname));
end
