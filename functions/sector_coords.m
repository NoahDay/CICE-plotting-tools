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
        coords = [-45,60;-75,60;-45,120;-75,120]; %(NW;SW,NE,SE) %% Need to verify these points with papers!!
    elseif sector == "WS" % Weddell sea
        coords = [];
    elseif sector == "SH" % Southern hemisphere
        coords = [-45,0;-89,0;-45,365;-89,365]; %(NW;SW,NE,SE)
    end
end
%------------- END OF CODE --------------