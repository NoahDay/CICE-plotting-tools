%sector_analysis analyse an Antarctic sector
%   There are three sectors available:
%       SA := South Africa, 
%       WS := Weddell Sea,
%       and EA := East Antarctica
% 
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
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon)

%% Ice Area
[lat,lon,row] = grid_read(grid);
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
data = data_format(filedir,'aice',row,lat,lon);

% Make sector mask
[len,wid] = size(data);
ocean_mask = data_format(filedir,'tmask',row,lat,lon);
sector_data = zeros(len,wid);
sector_mask = false(len,wid);
sector_data = ~ocean_mask*NaN;

for i = 0:lat_out(1)-lat_out(2)
    sector_mask(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = true;
    sector_data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3));
end


%data = data(sector_mask);
color_map = seaicecolormap();
latitude = [-90,-60];
longitude = [10,50];
figure(1)
w = worldmap('world');
    axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [0 0 0]);
    setm(w, 'maplatlimit', [-70,-40]);
    setm(w, 'maplonlimit', [10,50]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'on')
    setm(w, 'mlabellocation', 10);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', 0);
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,sector_data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %caxis([0,600]);
    %title("Representative FSD per cell (m) on 31-08-2005",'interpreter','latex','FontSize', 18)

%% Pancake ice concentration











    
    
%% Functions
function coords = sector_coords(sector)
% Coordinates of sector
%   There are three sectors available:
%       SA := South Africa, 
%       WS := Weddell Sea,
%       and EA := East Antarctica
    if sector == "SA"
        coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;NE,SW,SE)
    elseif sector == "EA"
        coords = [];
    elseif sector == "WS"
        coords = [];
    end
end