function [lat_out,lon_out] = lat_lon_finder(latit,longi,lat,lon)
% latit := desired latitude
% longi := desired longitude
% lat := latitude grid
% lon := longitude grid
%latit = -64

    idx = lat > latit;
    val = min(min(idx.*lat));
    idx = lon < longi;
    val2 = max(max(idx.*lon));

    pos = find(lat == val);
    pos2 = find(lon == val2);

    [len,wid] = size(lat);
    dummy_lat = zeros(len,wid);
    [len,wid] = size(lon);
    dummy_lon = zeros(len,wid);
    
    dummy_lat(pos) = lat(pos);
    dummy_lon(pos2) = lon(pos2);

    s = size(lon);
    [I,J] = ind2sub(s,find(dummy_lon));

    lon_out = I;%lon(I,J); 
    lat_out = find(dummy_lat(1,:));
end

