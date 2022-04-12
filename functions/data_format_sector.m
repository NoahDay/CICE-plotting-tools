function [sector_data, sector_mask] = data_format_sector(filedir,variable,sector)
%  if ~exist('user', 'var')
%         cd '/Users/noahday/GitHub/CICE-plotting-tools'
%  end
%  
%  if user == "a1724548"
%      cd '/Users/a1724548/GitHub/CICE-plotting-tools'
%  end
%cd '..'
addpath functions
    % Get the info
    info = ncinfo(filedir,variable);
    attributes = info.Attributes;
    coord_att = attributes(3); % Extract coordinate info
    coord_string = coord_att.Value;
    if coord_string(1:9) == 'TLON TLAT'
        coord_type = "t"; % T-grid
        coord_type = "t"; % T-grid
        lat = ncread(filedir,"TLAT");
        lon = ncread(filedir,"TLON");
    elseif coord_string(1:9) == 'ULON ULAT'
        coord_type = "u"; % U-grid
        lat = ncread(filedir,"ULAT");
        lon = ncread(filedir,"ULON");
    else
        error('ND: Unspecified grid type');
    end
    data_size = info.Size;
    idx = (data_size == 1);
    dim = length(data_size(~idx));
    if data_size(1) == 320 && data_size(2) == 384
        grid = "gx1";
        row = 37;
        lat = rearrange_matrix(lat,row,2);
        lon = rearrange_matrix(lon,row,2);

        lon = [zeros(1,384);lon];
        lat = [lat(1,:); lat];
    else
        error('ND: Unrecognised grid size.')
    end

    % Read in data
    data = data_format(filedir,variable);

    ocean_mask = data_format(filedir,'tmask');
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
                sector_data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3),:) = data(lon_out(1):lon_out(3), lat_out(1)-i:lat_out(3));
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
    else 
        error('ND: dimension not found')
    end
end