function [lat_ice_edge, lon_ice_edge] = find_ice_edge(aice,SIC,sector,lat,lon)
%FIND_ICE_EDGE finds the ice edge of ice area data aice according to a
%threshold SIC
%   Detailed explanation goes here
[len,wid] = size(aice);
ice_edge = zeros(len,wid);

coords = sector_coords(sector);
[lat_north, lon_east] = lat_lon_finder(coords(1,1),coords(1,2),lat,lon);
[lat_south, lon_west] = lat_lon_finder(coords(2,1),coords(3,2),lat,lon);

edge = ones(1,len);

for j = lon_east:lon_west
    for i = lat_north:-1:lat_south
        if aice(j,i) < SIC
            edge(j) = i+1;
        end
    end
end

for j = 1:len
    lat_ice_edge(j) = lat(j,edge(j));
    lon_ice_edge(j) = lon(j,edge(j));
end
end