filename = 'iced.2005-01-02-00000.nc';

% Read the header
ncdisp(filename)

u_vel=ncread(filename,'uvel');
v_vel=ncread(filename,'vvel');

% grid
lat = ncread('grid_gx3.nc','ulat');
lon = ncread('grid_gx3.nc','ulon');

% converting to degrees
deglon = rad2deg(lon);
deglat = rad2deg(lat);

[~, row] = patch_mat(deglon);
row = 90;
deglon = remove_row(deglon,row);
deglat = remove_row(deglat, row);


color_map = seaicecolormap();
%% Mapping

% sea ice area
sea_ice_area = ncread(filename, 'aicen');
si_area_01 = sea_ice_area(:,:,1);
si_area_02 = sea_ice_area(:,:,2);

latitude = [-90,90];
longitude = [-180,180];


%idx = si_area_01 == 0;
%si_area_01 = si_area_01 + 0*idx;

si_area_01 = remove_row(si_area_01, row);

w = worldmap(latitude,longitude);
axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
setm(w, 'Origin', [-90 0 0]);
setm(w, 'maplatlimit', [-90,-55]);
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


load coastlines
whos
[latcells, loncells] = polysplit(coastlat, coastlon);
numel(latcells)



%plotm(coastlat, coastlon)
%bedmap2('gl','xy')
pcolorm(deglat,deglon,si_area_01)
patchm(antarctica.Lat, antarctica.Lon, [1 1 1])
