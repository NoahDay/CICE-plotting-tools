function coords = sector_coords(sector)
% Coordinates of sector
%   There are three sectors available:
%       SA := South Africa, 
%       WS := Weddell Sea,
%       and EA := East Antarctica
    if sector == "SA"
        coords = [-55,20;-70,20;-55,40;-70,40]; %(NW;SW,NE,SE)
    elseif sector == "EA"
        coords = [-45,60;-75,60;-45,120;-75,120]; %(NW;SW,NE,SE) %% Need to verify these points with papers!!
    elseif sector == "WS"
        coords = [];
    elseif sector == "SH"
        coords = [-45,0;-89,0;-45,365;-89,365]; %(NW;SW,NE,SE)
    end
end