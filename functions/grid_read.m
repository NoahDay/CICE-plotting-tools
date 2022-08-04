 function [lat,lon,row] = grid_read(grid)
 dim = 2;
    if grid == "gx3"
        row = 11;

        ulat = ncread('grid/gx3/grid_gx3.nc','ulat');
        ulon = ncread('grid/gx3/grid_gx3.nc','ulon');

        % converting to degrees
        lon = rad2deg(ulon);
        lat = rad2deg(ulat);

        lat = rearrange_matrix(lat,row,dim);
        lon = rearrange_matrix(lon,row,dim);

    elseif grid == "om2"
    % Grid type from COSIMA
        lat = ncread('grid/om2/icegrid.nc','tlat')*180/pi;
        lon = mod(ncread('grid/om2/icegrid.nc','tlon')*180/pi,360);
        row = 1;
    else
        row = 37;
        lat = ncread('grid/gx1/global_gx1.bathy.nc','TLAT');
        lon = ncread('grid/gx1/global_gx1.bathy.nc','TLON');

        lat = rearrange_matrix(lat,row,dim);
        lon = rearrange_matrix(lon,row,dim);


        lon = [zeros(1,384);lon];
        lat = [lat(1,:); lat];
        
end
 end