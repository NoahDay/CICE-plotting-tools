function [lat_out,lon_out] = lat_lon_finder(latit,longi,lat,lon)
% latit := desired latitude (-90,90)
% longi := desired longitude (0,360)
% lat := latitude grid
% lon := longitude grid
%latit = -64
    
    val = min(min(abs(lat - latit)));%min(min(idx.*lat));
    pos = find(abs(lat - latit) == val);
    [len,wid] = size(lat);
    dummy_lat = zeros(len,wid);
    dummy_lat(pos) = lat(pos);
    lat_out = find(dummy_lat(1,:));
    lon_vec = lon(:,lat_out);
    val2 = min(abs(lon_vec - longi));
    lon_out = find(abs(lon_vec - longi) == val2);     
end

