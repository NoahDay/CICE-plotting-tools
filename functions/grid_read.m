function [lat,lon,row,ulat,ulon] = grid_read(grid)
    cd '/Users/noahday/GitHub/CICE-plotting-tools/'
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
        tlat = ncread('grid/om2/icegrid.nc','tlat')*180/pi;
        tlon = mod(ncread('grid/om2/icegrid.nc','tlon')*180/pi,360);
        ulat = ncread('grid/om2/icegrid.nc','ulat')*180/pi;
        ulon = mod(ncread('grid/om2/icegrid.nc','ulon')*180/pi,360);
        row = 1;
        %elseif data_size(1) == 360 && data_size(2) == 300
        row = 281;
        % OM2 grid
        tlat = rearrange_matrix(tlat,row,2);
        tlon = rearrange_matrix(tlon,row,2);
        lon = [zeros(1,300);tlon];
        lat = [tlat(1,:); tlat];

        row = 280;
        ulat = rearrange_matrix(ulat,row,2);
        ulon = mod(rearrange_matrix(ulon,row,2)+360,360);
        ulon = [zeros(1,300);ulon];
        ulat = [ulat(1,:); ulat];
    else
        row = 37;
        lat = ncread('grid/gx1/global_gx1.bathy.nc','TLAT');
        lon = ncread('grid/gx1/global_gx1.bathy.nc','TLON');

        lat = rearrange_matrix(lat,row,dim);
        lon = rearrange_matrix(lon,row,dim);


        lon = [zeros(1,384);lon];
        lat = [lat(1,:); lat];
        

        ulat = ncread('grid/gx1/global_gx1.bathy.nc','ULAT');
        ulon = ncread('grid/gx1/global_gx1.bathy.nc','ULON');

        ulat = rearrange_matrix(ulat,row,dim);
        ulon = rearrange_matrix(ulon,row,dim);


        ulon = [zeros(1,384);ulon];
        ulat = [ulat(1,:); ulat];
end
 end