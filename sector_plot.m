%sector_plot plot the chosen Antarctic sector.
%   There are two sectors available, SA := South Africa, and EA := East
%   Antarctic
%% Read in the data.
clear
close all
addpath functions
filename = 'cases/init/history/iceh_ic.2005-01-01-03600.nc';
%ncdisp(filename)

% Parameters
sector = "SA";
grid = 'gx1';
filedir = "cases/init/history/iceh_ic.2005-01-01-03600.nc";

% Coordinates
if sector == "SA"
    coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
elseif sector == "EA"
    coords = [];
%% Retreive the grid.
[lat,lon,row] = grid_read(grid);
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
data = data_format(filedir,'aice',row,lat,lon);
% Vertices
for i = 1:4
    data(lon_out(i),lat_out(i)) = 10;
end
% Edges
% North latitude
data(lon_out(1):lon_out(3), lat_out(1):lat_out(3)) = 10;
% South latitude
data(lon_out(2):lon_out(4), lat_out(2):lat_out(4)) = 10;
% East longitude
data(lon_out(2):lon_out(1), lat_out(2):lat_out(1)) = 10;
% West longitude
data(lon_out(4):lon_out(3), lat_out(4):lat_out(3)) = 10;
%% Begin mapping.
color_map = seaicecolormap();
latitude = [-90,-60];
longitude = [10,50];
figure(1)
w = worldmap('world');
    axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [0 0 0]);
    setm(w, 'maplatlimit', [-90,-40]);
    setm(w, 'maplonlimit', [10,50]);
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

