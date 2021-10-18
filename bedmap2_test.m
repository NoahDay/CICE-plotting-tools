close all
shelf = bedmap2_data('icemask');
shelf = 10.0*(shelf == 1); 
shelf(shelf==0) = NaN;% 1 iceshelf
[latshelf,lonshelf] = bedmap2_data('latlon');
%pcolorps(lat,lon,shelf)

filename = 'casew/history/iceh.2005-01-10.nc';

% Read the header
ncdisp(filename);
% 
% u_vel=ncread(filename,'uvel');
% v_vel=ncread(filename,'vvel');

% grid
lat = ncread('global_gx1.bathy.nc','TLAT');
lon = ncread('global_gx1.bathy.nc','TLON');

% converting to degrees
%deglon = rad2deg(lon);
%deglat = rad2deg(lat);

%row = 80;
%copies = 5;

lat = rearrange_matrix(lat,37);
lon = rearrange_matrix(lon,37);


lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

% sea ice area
data = ncread(filename, 'aice');
%data = data(:,:,2);
data = rearrange_matrix(data,37);
data = [data; data(end,:)];


color_map = seaicecolormap();
%% Mapping



latitude = [-90,90];
longitude = [-180,180];


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
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    caxis([0 1])
    colorbar
    pcolorm(latshelf,lonshelf,shelf)
    
    