clear all
close all
addpath functions
filename = 'DATA/gx1v3/forcing/prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_195801010130-195812312230.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
% Read the header
ncdisp(filename)
% 
grid = "gx1";
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
%[lat,lon,row] = grid_read(grid);
data = ncread(filename, 'prra');
data = data(:,:,1);
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