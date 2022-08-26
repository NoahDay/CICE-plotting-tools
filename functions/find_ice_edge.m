function [lat_ice_edge_new, lon_ice_edge_new, edge] = find_ice_edge(aice,SIC,sector,lat,lon)
%FIND_ICE_EDGE finds the ice edge of ice area data aice according to a
%threshold SIC
%   Detailed explanation goes here
    [len,wid] = size(aice);
    ice_edge = zeros(len,wid);
    
    coords = sector_coords(sector);
    [lat_north, lon_east] = lat_lon_finder(coords(1,1),coords(1,2),lat,lon);
    [lat_south, lon_west] = lat_lon_finder(-75,coords(3,2),lat,lon);

    edge = ones(1,len);
    
    if lon_east > lon_west
        for j = lon_east:len
            for i = lat_north:-1:lat_south
                if aice(j,i) > SIC
                    edge(j) = i+1;
                    break
                end
            end
        end
        for j = 1:lon_west
            for i = lat_north:-1:lat_south
                if aice(j,i) > SIC
                    edge(j) = i+1;
                    break
                end
            end
        end

    else
        for j = lon_east:lon_west
            for i = lat_north:-1:lat_south
                if aice(j,i) > SIC
                    edge(j) = i+1;
                    break
                end
            end
        end
    end
    
    for j = 1:len
        lat_ice_edge(j) = lat(j,edge(j));
        lon_ice_edge(j) = lon(j,edge(j));
    end
    
    % Make the ice edge fit the resolution
    count = 2;
    lat_ice_edge_new(1) = lat_ice_edge(1);
    lat_ice_edge_new(1) = lat_ice_edge(1);
    for j = 2:len
        if lat_ice_edge(j) ~= lat_ice_edge(j-1)
            % Then we have a step
            lat_ice_edge_new(count) = lat_ice_edge(j);
            lon_ice_edge_new(count) = lon_ice_edge(j-1);
            count = count + 1;
        end
        lat_ice_edge_new(count) = lat_ice_edge(j);
        lon_ice_edge_new(count) = lon_ice_edge(j);
        count = count + 1;
    end

end