clear all
close all
filedir =  'grid/ww3_200506.nc'%'grid/gridded_ww3.glob_24m.200501.nc';%'ww3_gx1.nc';%'gridded_ww3.glob_24m.200501.nc';%'ww3.20100101_efreq_remapgx3.nc';%'gridded_ww3.glob_24m.200501.nc';%'ww3.200501_spec.nc';

ncdisp(filedir)

% grid
%lat = ncread('global_gx1.bathy.nc','TLAT');
%lon = ncread('global_gx1.bathy.nc','TLON');
lat = ncread('ww3.200501_spec.nc', 'latitude');
lon = ncread('ww3.200501_spec.nc', 'longitude');


% converting to degrees
%deglon = rad2deg(lon);
%deglat = rad2deg(lat);

%row = 80;
%copies = 5;

 row = 11;
    
    %ulat = ncread('grid_gx3.nc','ulat');
    %ulon = ncread('grid_gx3.nc','ulon');

    % converting to degrees
 %   lon = rad2deg(ulon);
  %  lat = rad2deg(ulat);

   % lat = rearrange_matrix(lat,row);
    %lon = rearrange_matrix(lon,row);
   

% sea ice area
data = ncread(filedir, 'hs');% 'efreq'); % wave_sig_ht, dafsd_wave, fsdrad, peak_period, aice
size(data)
data = data(:,:,500);%*pi/180;
%data(:,48) = 10;
%data(1:10,48) = 5;
%data = rearrange_matrix(data,37);
%data = [data; data(end,:)];

%lat =  ncread(filedir, 'latitude');
%lon = ncread(filedir, 'longitude');

%% Mapping
color_map = seaicecolormap();


latitude = [-78,78];
longitude = [0,360];


w = worldmap('world');
axesm eqdazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
setm(w, 'Origin', [-90 0 0]);
setm(w, 'maplatlimit', [-90,-40]);
%setm(w, 'maplonlimit', [0,360]);
%setm(w, 'meridianlabel', 'on')
%setm(w, 'parallellabel', 'off')
%setm(w, 'mlabellocation', 30);
%setm(w, 'plabellocation', 10);
%setm(w, 'mlabelparallel', -55);
%setm(w, 'grid', 'on');
%setm(w, 'labelrotation', 'on')

pcolorm(latitude,longitude,data')
colorbar
%caxis([0,10]);
title("Mean wave direction")


