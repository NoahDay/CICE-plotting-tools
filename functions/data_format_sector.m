function [sector_data, sector_mask] = data_format_sector(filedir,variable,sector,dim)
%  if ~exist('user', 'var')
%         cd '/Users/noahday/GitHub/CICE-plotting-tools'
%  end
%  
%  if user == "a1724548"
%      cd '/Users/a1724548/GitHub/CICE-plotting-tools'
%  end
%cd '..'
addpath functions
     if ~exist('dim', 'var')
        dim = 2; 
     end

    % Grid
    grid = "gx1";
    [lat,lon,row] = grid_read(grid);
    % Read in data
    data = data_format(filedir,variable,row,lat,lon,dim);
    ocean_mask = data_format(filedir,'tmask',row,lat,lon);
    if dim == 2
        if isstring(sector)
            % Export a grid of the the specific sector
            % Coordinates
            coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
            for i = 1:4
                [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
            end
            [len,wid] = size(data);
            sector_data = zeros(len,wid);
            sector_mask = false(len,wid);
            sector_data = ~ocean_mask*NaN;
            for i = 0:lat_out(1)-lat_out(2)
                sector_mask(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = true;
                sector_data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3));
            end
        else % Transect
            [len,wid] = size(data);
            [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            sector_mask = false(len,wid);
            sector_mask(lon_out,:) = true;
            sector_data(lon_out,:) = data(lon_out,:);
        end
    elseif dim == 3
        if isstring(sector)
            % Export a grid of the the specific sector
            % Coordinates
            coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
            for i = 1:4
                [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
            end
            [len,wid,dep] = size(data);
            sector_data = zeros(len,wid,dep);
            sector_mask = false(len,wid);
            %sector_data = ~ocean_mask*NaN;
            for i = 0:lat_out(1)-lat_out(2)
                sector_mask(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = true;
                sector_data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3),:) = data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3),:);
            end
        else % Transect
            [len,wid] = size(data);
            [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            sector_mask = false(len,wid);
            sector_mask(lon_out,:) = true;
            sector_data(lon_out,:,:) = data(lon_out,:,:);
        end
    elseif dim == 4
        if isstring(sector)
            % Export a grid of the the specific sector
            % Coordinates
            coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
            for i = 1:4
                [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
            end
            [len,wid,dep,hei] = size(data);
            sector_data = zeros(len,wid,dep,hei);
            sector_mask = false(len,wid);
            %sector_data = ~ocean_mask*NaN;
            for i = 0:lat_out(1)-lat_out(2)
                sector_mask(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3)) = true;
                sector_data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3),:,:) = data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3),:,:);
            end
        else % Transect
            [len,wid] = size(data);
            [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            sector_mask = false(len,wid);
            sector_mask(lon_out,:) = true;
            sector_data(lon_out,:,:,:) = data(lon_out,:,:,:);
        end
    end
end