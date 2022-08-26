close all
addpath packages/bedmap2_toolbox_v4

shelf = bedmap2_data('icemask');
shelf = 0.0*(shelf == 1); 
shelf(shelf==0) = NaN; % 1 iceshelf
[latshelf,lonshelf] = bedmap2_data('latlon');

%%
filename = '/Users/noahday/GitHub/CICE-plotting-tools/bedmap2shelf.mat';
 Z = double(imread(filename));

%%
filename = 'cases/8month/history/iceh.2005-01-10.nc';

% Read the header
ncdisp(filename);

% grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');
lat = rearrange_matrix(lat,37,2);
lon = rearrange_matrix(lon,37,2);
lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

% sea ice area
data = ncread(filename, 'aice');
%data = data(:,:,2);
data = rearrange_matrix(data,37,2);
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
    