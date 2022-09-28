function coords = sector_coords(sector)
%SECTOR_COORDS - Retrieves the coordinates of predetermined sectors.
%
% Syntax:  [coordinates] = sector_coords(sector)
%
% Inputs:
%    sector - string
%
% Outputs:
%    coordinates - 4x2 matrix
%
% Example: 
%    coords = sector_coords("SH")
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Noah Day
% email: noah.day@adelaide.edu.au
% Mar 2022; Last revision: 25-March-2022

%------------- BEGIN CODE --------------
    if sector == "SA" % South African sector
        coords = [-55,20;-70,20;-55,40;-70,40]; %(NW;SW,NE,SE)
    elseif sector == "EA" % Eastern Antarctica sector
        coords = [-45,1;-75,1;-45,140;-75,140]; %(NW;SW,NE,SE) %% Need to verify these points with papers!!
    elseif sector == "AU" % Australia
        coords = [-25,1;-75,1;-25,150;-75,150]; %(NW;SW,NE,SE) %% Need to verify these points with papers!!
    elseif sector == "WS" % Weddell sea
        coords = [];
    elseif sector == "SH" % Southern hemisphere
        coords = [-25,0;-89,0;-25,365;-89,365]; %(NW;SW,NE,SE)
    elseif sector == "vichi" % Coords from Vichi et al. (2019) Figure 1.
        max_lat = -40;
        min_lon = 365-20;
        max_lon = 40;
        min_lat = -75;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    elseif sector == "vichi1" % Coords from Vichi et al. (2019) Figure 1.
        max_lat = -40;
        min_lon = 365-20;
        %max_lon = 40;
        max_lon = 365;
        min_lat = -65;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    elseif sector == "vichi2" % Coords from Vichi et al. (2019) Figure 1.
        max_lat = -40;
        min_lon = 1;
        max_lon = 60;
        min_lat = -75;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    elseif sector == "vichi3" % Coords from Vichi et al. (2019) Figure 3.
        max_lat = -55;
        min_lon = 5;
        max_lon = 40;
        min_lat = -75;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    elseif sector == "toyota2016" % Coords from Vichi et al. (2019) Figure 3.
        max_lat = -55;
        min_lon = 105;
        max_lon = 130;
        min_lat = -75;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    elseif sector == "global" % Coords from Vichi et al. (2019) Figure 3.
        max_lat = 0;
        min_lon = 0;
        max_lon = 365;
        min_lat = -90;
        coords = [max_lat,min_lon;min_lat,min_lon;max_lat,max_lon;min_lat,max_lon]; %(NW;SW,NE,SE)
    else
        error("ND: no sector specified")
    end
end
%------------- END OF CODE --------------