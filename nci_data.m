% Read and plot NCI JRA55 data

clear
close all
addpath functions

% Parameters
%sector = "EA";
grid = "om2_1deg";
[lat,lon,row] = grid_read(grid);
filedir = '/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/om2_1deg/JRA55-do-1-5-0/atmos/3hr/prra/prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_200501010130-200512312230.nc';

ncdisp(filedir)
lat = ncread(filedir,'lat');
lon = ncread(filedir,'lon');
prra = ncread(filedir,'prra');
prra_1 = prra(:,:,1);

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
    pcolorm(lat,lon,prra_1)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %%
fileID = fopen('grid_gx1.bin');
A = fread(fileID);
fclose(fileID);