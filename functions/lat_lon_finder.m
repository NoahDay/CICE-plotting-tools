function [lat_out,lon_out] = lat_lon_finder(latit,longi,lat,lon)
% latit := desired latitude
% longi := desired longitude
% lat := latitude grid
% lon := longitude grid
%latit = -64
    idx = lat > latit;
    val = min(min(idx.*lat));
    pos = find(lat == val);
    [len,wid] = size(lat);
    dummy_lat = zeros(len,wid);
    dummy_lat(pos) = lat(pos);
    lat_out = find(dummy_lat(1,:));
    lon_vec = lon(:,lat_out);
    idx = lon_vec < longi;
    val2 = max(max(idx.*lon_vec));
    lon_out = find(lon_vec == val2);     
end

