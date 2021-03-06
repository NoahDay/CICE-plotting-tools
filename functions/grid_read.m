 function [lat,lon,row] = grid_read(grid)
 dim = 2;
    if grid == 'gx3'
        row = 11;

        ulat = ncread('grid/gx3/grid_gx3.nc','ulat');
        ulon = ncread('grid/gx3/grid_gx3.nc','ulon');

        % converting to degrees
        lon = rad2deg(ulon);
        lat = rad2deg(ulat);

        lat = rearrange_matrix(lat,row,dim);
        lon = rearrange_matrix(lon,row,dim);

    elseif grid == "om2_1deg"
    % Grid type from COSIMA
        lat = ncread('grid/om2_1deg/om2_1deg_basinmask_yt_xt.nc','yt_ocean');
        lon = ncread('grid/om2_1deg/om2_1deg_basinmask_yt_xt.nc','xt_ocean');
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