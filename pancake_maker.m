function y = pancake_maker(filename1, plot_title, variable, grid, map_type, user, rad)
% Plot the pancake ice

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

% data 1
data = ncread(filename1, variable);
[~, ~, n] = size(data);

for level = 1:n
    data_1 = data(:,:,level);
end


    latitude = [-90,90];
    longitude = [-180,180];

    data_1 = rearrange_matrix(data_1,row);


    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_1 = [data_1; data_1(end,:)];
    
    % setting the maximum radius 
    
    idx = data_1 < rad;
    data = data_1.*idx;

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
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

    %antarctica = shaperead('landareas', 'UseGeoCoords', true,...
    %  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

  
    
    title(plot_title);
    %caxis([-10 10]) %fsdrad [80 250]
    colorbar
    
end