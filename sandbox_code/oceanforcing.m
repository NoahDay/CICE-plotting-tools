clear all
close all
addpath functions
%cd .. % Move up one directory
filename = 'ocean_forcing_clim_2D_gx1.20210330.nc';
% '/Users/noahday/Gadi/2010/JRA55_03hr_forcing_2010.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
%filename = 'DATA/CESM/MONTHLY/ocean_forcing_clim_2D_gx1.20210330.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';
%lat = ncread('grid/gx1/global_gx1.bathy.nc','TLAT');
%lon = ncread('grid/gx1/global_gx1.bathy.nc','TLON');
%[lat,lon,row] = grid_read(grid);
%data = ncread(filename, 'U');
%data = data_format(filename,'prra',row,lat,lon,3);
%data = data(:,:,1);
[lat,lon,row] = grid_read(grid);
lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);
dim = 3;

data_u = data_format(filename,'U',row,lat,lon,dim);
data_u2 = data_u(:,:,1);
data_v = data_format(filename,'V',row,lat,lon,dim);
data_v2 = data_v(:,:,1);
u_vec = reshape(data_u2,1,[]);
v_vec = reshape(data_v2,1,[]);

data2 = data_format(filename,'hblt',row,lat,lon,dim);%sqrt(data_u2.^2 + data_v2.^2);
data = data2(:,:,1);
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    caxis([0,100])

%%
clear all; close all; clc
addpath functions
filename = "/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh.2005-01-01.nc";
[lat_cice,lon_cice] = grid_read("om2");

sst = data_format(filename, "vocn");
conFigure(30)
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat_cice,lon_cice,sst)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    caxis([0,1])


%%
close all; clc
addpath functions
filename = "/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/ocean/output2010.nc";
%[lat,lon] = grid_read("om2");
tmask = ncread("/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh.2005-01-01.nc",'tmask');
data = ncread(filename, "u");
datav = ncread(filename, "v");
lat = ncread(filename, "yu_ocean");
lon = ncread(filename, "xu_ocean");
data1d(:,:) = data(:,:,1);
datav1d(:,:) = datav(:,:,1);

conFigure(30)
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat',lon',(data1d.*tmask)')
    land = shaperead('landareas', 'UseGeoCoords', true);
    %geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    %caxis([120,130])

    min(min(data1d))

sprintf('The maximum u-velocity is: %g', max(max(data1d)))
sprintf('The minimum v-velocity is: %g', min(min(datav1d)))

sprintf('The maximum v-velocity is: %g', max(max(datav1d)))
sprintf('The minimum v-velocity is: %g', min(min(datav1d)))



%% JRA55 data
close all; clc
addpath functions
filename = "/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/ocean/output2010.nc";
%[lat,lon] = grid_read("om2");
tmask = ncread("/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh.2005-01-01.nc",'tmask');
data = ncread(filename, "u");
datav = ncread(filename, "v");
lat = ncread(filename, "yu_ocean");
lon = ncread(filename, "xu_ocean");
data1d(:,:) = data(:,:,1);
datav1d(:,:) = datav(:,:,1);

conFigure(30)
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
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat',lon',(data1d.*tmask)')
    land = shaperead('landareas', 'UseGeoCoords', true);
    %geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    %caxis([120,130])

    min(min(data1d))

sprintf('The maximum u-velocity is: %g', max(max(data1d)))
sprintf('The minimum v-velocity is: %g', min(min(datav1d)))

sprintf('The maximum v-velocity is: %g', max(max(datav1d)))
sprintf('The minimum v-velocity is: %g', min(min(datav1d)))



