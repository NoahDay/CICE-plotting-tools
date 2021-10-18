clear all
close all
filename = 'cases/20day/history/iceh.2005-01-18.nc';

% Read the header
ncdisp(filename)
% 
% u_vel=ncread(filename,'uvel');
% v_vel=ncread(filename,'vvel');

% grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');

% converting to degrees
%deglon = rad2deg(lon);
%deglat = rad2deg(lat);

%row = 80;
%copies = 5;
dim = 2;
lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);


lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];



% sea ice area
data = ncread(filename, 'wave_sig_ht'); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
%data = data(:,:,6);
data = rearrange_matrix(data,37,dim);
data = [data; data(end,:)];


color_map = seaicecolormap();
%% Mapping



latitude = [-90,-30];
longitude = [-180,180];


%idx = si_area_01 == 0;
%si_area_01 = si_area_01 + 0*idx;

%si_area_01 = rearrange_matrix(si_area_01,11);


%deglon = [deglon; deglon(100,:)+3.6];
%deglat = [deglat; deglat(100,:)];
%si_area_01 = [si_area_01; si_area_01(100,:)];


w = worldmap('world');
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
%land = shaperead('landareas', 'UseGeoCoords', true);
%geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

%antarctica = shaperead('landareas', 'UseGeoCoords', true,...
  %'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});


%load coastlines
%whos
%[latcells, loncells] = polysplit(coastlat, coastlon);
%numel(latcells)


%lon = linspace(0,360,100);
%lon = repmat(lon,116,1)';
%plotm(coastlat, coastlon)
%bedmap2('gl','xy')
%bedmap2('shelves')
pcolorm(lat,lon,data)
%quiverm(lat,lon,data)
lat0 = [-90 -39.7]; 
lon0 = [-5.4 2.9];
deltalat = [5 5]; 
deltalon = [3 3];
%quiverm(lat0,lon0,deltalat,deltalon,'r')
colorbar
%caxis([0,10]);
%title("Mean wave direction")
%caxis([-1 1])
%patchm(antarctica.Lat, antarctica.Lon, [1 1 1])
%saveas(gcf,'test.png')