%Initial_ice is a report of the ice initial conditions for CICE. 
%   It is taken from point observations in Perovich et al. (2014) from the
%   Arctic. This investigation will be targetted at the FSD, ITD, ice 
%   extent around Antarctica.
%% Read in the data.
clear
close all
addpath functions
filedir = 'cases/init/history/iceh_ic.2005-01-01-03600.nc';
% Files
% Initial conditions forcing: iced_gx1_v6.2005-01-01.nc
% Initial conditions of CICE: iceh_ic.2005-01-01-03600.nc
% First timestep (1 hour): iceh_inst.2005-01-01-03600.nc
% After one day: iceh.2005-01-01.nc
ncdisp(filedir)

% Parameters
fstd_switch = 0;
coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
grid = 'gx1';
%% Retreive the grid.
[lat,lon,row] = grid_read(grid);
% for i = 1:4
%     [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
% end
if fstd_switch == 1
    data = data_format(filedir,'aicen_d',row,lat,lon,4);
    [lat_out,lon_out] = lat_lon_finder(coords(2,1),coords(2,2),lat,lon);
    fstd = data(lon_out,lat_out,:,:);
    size(sum(fstd))
end

%% Begin mapping.
figure(1)
data = data_format(filedir,'afsd',row,lat,lon,3); %aicen
for i = 1:24
    idx = data(:,:,i)> 0;
    sum(sum(idx))
end
%data(lon_out,lat_out) = 10;
data = data(:,:,10);
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

