clear all
close all
addpath functions
day = 6;
month = 1;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end
case_name = 'ocntest';
ticker = 1;
SIC = 0.15; 
filename = strcat('cases/',case_name,"/history/iceh.",date,".nc");
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
dim = 2;

data_u = data_format(filename,'strocnx',row,lat,lon,dim);
data_u2 = data_u(:,:,1);
data_v = data_format(filename,'strocny',row,lat,lon,dim);
data_v2 = data_v(:,:,1);
u_vec = reshape(data_u2,1,[]);
v_vec = reshape(data_v2,1,[]);

data = sqrt(data_u2.^2 + data_v2.^2);
data =  data_format(filename,'dvidtd',row,lat,lon,dim);
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