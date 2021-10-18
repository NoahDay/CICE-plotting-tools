function y = map_creator(filename, plot_title_vec, i, variable, grid, map_type, user)
%MAP_CREATOR creates map images of Antarctica given netcdf data files
%   filename: the directory to the .nc file
%   plot_title: string containing the title for the plot
%   i: the date of the data file
%   variable: the variable in the dataset we want to plot
%   grid: specify the grid. eg. 'gx1', 'gx3'

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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_1)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

    %antarctica = shaperead('landareas', 'UseGeoCoords', true,...
    %  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

 if variable == "fsdrad"
         plot_variable = "FSD radius ";
     elseif variable == "wave_sig_ht"
         plot_variable = "Significant wave height ";
     elseif variable == "peak_period"
         plot_variable = "Peak period ";
     elseif variable == "aice"
         plot_variable = "Concentration of ice ";
     elseif variable == "mean_wave_dir"
         plot_variable = "Mean wave direction (rads) ";
     else
         plot_variable = variable;
 end
    fontSize = 14; 
    plot_title = strcat(plot_variable, plot_title_vec(i));
    title(plot_title, 'FontSize', fontSize);
    %caxis([0 400]) %fsdrad [80 250]
    colorbar
    limit = colorlims(variable);
    caxis(limit);
    %patchm(antarctica.Lat, antarctica.Lon, [1 1 1]);
    figname = sprintf('image%d.png', i); 
    text = sprintf('/Users/%s/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames', user);
    fname = '/Volumes/Noah SSD/MATLAB/PhD Project/CICE Plotting/frames'; 
    %'/Users/noahday/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames';
    saveas(gcf,fullfile(fname, figname));
end
end

