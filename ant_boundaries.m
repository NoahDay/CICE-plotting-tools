clear all
close all
filename = 'cases/4proc/history/iceh.2005-01-01.nc';

% Read the header
ncdisp(filename)
dim = 2;

% grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');


lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);


lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];



% sea ice area
data = ncread(filename, 'wave_sig_ht'); % wave_sig_ht, dafsd_wave, fsdrad
%data = data(:,:,5);
data = rearrange_matrix(data,37,dim);
data = [data; data(end,:)];
color_map = seaicecolormap();
%% Mapping



latitude = [-90,90];
longitude = [-180,180];


% w = worldmap('world');
% axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
% setm(w, 'Origin', [-90 0 0]);
% setm(w, 'maplatlimit', [-90,-40]);
% setm(w, 'maplonlimit', [-180,180]);
% setm(w, 'meridianlabel', 'on')
% setm(w, 'parallellabel', 'off')
% setm(w, 'mlabellocation', 30);
% setm(w, 'plabellocation', 10);
% setm(w, 'mlabelparallel', -55);
% setm(w, 'grid', 'on');
% setm(w, 'labelrotation', 'on')

%land = shaperead('landareas', 'UseGeoCoords', true);
%geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])



% pcolorm(lat,lon,data)
% colorbar
% hold on
[lat2,lon2] = antbounds_data('coast');
plotps(lat2,lon2,'blue')
[gllat,gllon] = antbounds_data('gl');
plotps(gllat,gllon,'red')

% [lat2,lon2] = antbounds_data(datatype);
% [x,y] = antbounds_data(datatype,'xy')
% [...,names] = antbounds_data('shelves')
% plotps(lat2,lon2,'.','color',0.8*[1 1 1]);
% hold on
% axis tight
% %antbounds('gl','black')
% %antbounds('coast','black')
% %mapzoomps('nw')
% hold off



