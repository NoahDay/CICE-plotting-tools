%Initial_ice is a report of the ice initial conditions for CICE. 
%   It is taken from point observations in Perovich et al. (2014) from the
%   Arctic. This investigation will be targetted at the FSD, ITD, ice 
%   extent around Antarctica.
%% Read in the data.
clear
close all
addpath functions
filedir = 'cases/init/history/iceh_ic.2005-01-01-03600.nc';
ncdisp(filedir)

% Parameters
coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
grid = 'gx1';
%% Retreive the grid.
[lat,lon,row] = grid_read(grid);
% for i = 1:4
%     [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
% end
data = data_format(filedir,'aicen_d',row,lat,lon,4);
[lat_out,lon_out] = lat_lon_finder(coords(2,1),coords(2,2),lat,lon);
fstd = data(lon_out,lat_out,:,:);
size(sum(fstd))

% % Vertices
% for i = 1:4
%     data(lon_out(i),lat_out(i)) = 10;
% end
% % Edges
% % North latitude
% data(lon_out(1):lon_out(3), lat_out(1):lat_out(3)) = 10;
% % South latitude
% data(lon_out(2):lon_out(4), lat_out(2):lat_out(4)) = 10;
% % East longitude
% data(lon_out(2):lon_out(1), lat_out(2):lat_out(1)) = 10;
% % West longitude
% data(lon_out(4):lon_out(3), lat_out(4):lat_out(3)) = 10;
%% Begin mapping.
figure(1)
data = data_format(filedir,'fsdrad',row,lat,lon,4);
data(lon_out,lat_out) = 10;
color_map = seaicecolormap();
latitude = [-90,-60];
longitude = [10,50];
figure(1)
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-40]);
    %setm(w, 'maplonlimit', [10,50]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    %setm(w, 'mlabellocation', 30);
    %setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', 0);
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %caxis([0,600]);
    %title("Representative FSD per cell (m) on 31-08-2005",'interpreter','latex','FontSize', 18)

