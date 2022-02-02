function coords = sector_coords(sector)
% Coordinates of sector
%   There are three sectors available:
%       SA := South Africa, 
%       WS := Weddell Sea,
%       and EA := East Antarctica
    if sector == "SA"
        coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
    elseif sector == "EA"
        coords = [];
    elseif sector == "WS"
        coords = [];
    end
end