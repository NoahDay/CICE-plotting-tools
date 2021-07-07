function y = map_creator(filename, plot_title, i, variable, grid)
%MAP_CREATOR creates map images of Antarctica given netcdf data files
%   filename: the directory to the .nc file
%   plot_title: string containing the title for the plot
%   i: the date of the data file
%   variable: the variable in the dataset we want to plot
%   grid: specify the grid. eg. 'gx1', 'gx3'

% grid
if grid == 'gx3'
    row = 11;
    
    ulat = ncread('grid_gx3.nc','ulat');
    ulon = ncread('grid_gx3.nc','ulon');

    % converting to degrees
    lon = rad2deg(ulon);
    lat = rad2deg(ulat);

    lat = rearrange_matrix(lat,row);
    lon = rearrange_matrix(lon,row);
   
else
    row = 37;
    lat = ncread('global_gx1.bathy.nc','TLAT');
    lon = ncread('global_gx1.bathy.nc','TLON');

    lat = rearrange_matrix(lat,row);
    lon = rearrange_matrix(lon,row);


    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
        
end


data = ncread(filename, variable);
[~, ~, n] = size(data);

for level = 1:n
    data_1 = data(:,:,level);

    latitude = [-90,90];
    longitude = [-180,180];

    data_1 = rearrange_matrix(data_1,row);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_1 = [data_1; data_1(end,:)];

    %% Mapping
    color_map = seaicecolormap();

    w = worldmap(latitude,longitude);
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-40]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -55);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')

    antarctica = shaperead('landareas', 'UseGeoCoords', true,...
      'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

    pcolorm(lat,lon,data_1);
    title(plot_title);
    colorbar;
    patchm(antarctica.Lat, antarctica.Lon, [1 1 1]);
    figname = sprintf('image%d.png', i);
    fname = '/Users/noahday/MATLAB-Drive/MATLAB/PhD Project/frames';
    saveas(gcf,fullfile(fname, figname));
end
end

