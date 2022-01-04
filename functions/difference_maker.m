function y = difference_maker(filename1, filename2, plot_title, variable, grid, map_type, user)
%Data 1 - Data 2
dim = 2;
% grid
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

% data 1
v1 = strcat(variable);%,'_d');
data = ncread(filename1, v1);
[~, ~, n] = size(data);

for level = 1:n
    data_1 = data(:,:,level);
end
% data 2
data2 = ncread(filename2, variable);
[~, ~, n] = size(data2);

for level = 1:n
    data_2 = data2(:,:,level);
end


    latitude = [-90,90];
    longitude = [-180,180];

    data_1 = rearrange_matrix(data_1,row,dim);
    data_2 = rearrange_matrix(data_2,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_1 = [data_1; data_1(end,:)];
    data_2 = [data_2; data_2(end,:)];
    
    data_1 = data_1 - data_2;

    %% Mapping
    
    color_map = seaicecolormap();
    w = worldmap('world');
        if map_type == "cassini"
            x_origin = -20;
            axesm cassini;
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlabellocation', 30);
        elseif map_type == "eqdcylin"
            x_origin = 0;
            axesm eqdcylin;
            setm(w, 'mlabelparallel', -85);
            setm(w, 'mlabellocation', 60);
        else
            x_origin = -90;
            axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlabellocation', 30);
        end
        setm(w, 'Origin', [x_origin 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_1)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

    %antarctica = shaperead('landareas', 'UseGeoCoords', true,...
    %  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});

  
    
    title(plot_title);

    %caxis(limits)
    colorbar
    
end