 function data_out = data_format(filedir,variable)
 % Extract the data from a netCDF file.
 % Input: (filedir, variable)
 % Type: (string, string)
 %
 % filedir := file directory
 % variable := variable of data we want to extract
 % dim := Dimension of variable
 % lat := grid latitude
 % lon := grid longitude

info = ncinfo(filedir,variable);
attributes = info.Attributes;
coord_att = attributes(3); % Extract coordinate info
coord_string = coord_att.Value;
if coord_string(1:9) == 'TLON TLAT'
    coord_type = "t"; % T-grid
    lat = ncread(filedir,"TLAT");
    lon = ncread(filedir,"TLON");
else
    coord_type = "u"; % U-grid
    lat = ncread(filedir,"ULAT");
    lon = ncread(filedir,"ULON");
end
data_size = info.Size;
dim = length(data_size);
if data_size(1) == 320 && data_size(2) == 384
    grid = "gx1";
    row = 37;
    lat = rearrange_matrix(lat,row,2);
    lon = rearrange_matrix(lon,row,2);

    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
elseif data_size(1) == 360 && data_size(2) == 300
        % OM2 grid
         if coord_type == "u"
            row = 280;
            lat = rearrange_matrix(lat,row,2);
            lon = mod(rearrange_matrix(lon,row,2)+360,360);
        else % t-grid
           % OM2 grid
            row = 281;
            lat = rearrange_matrix(lat,row,2);
            lon = rearrange_matrix(lon,row,2);
        end
end

 
 if dim == 3
    data1 = ncread(filedir, variable);

    latitude = [-90,90];
    longitude = [-180,180];

    data1 = rearrange_matrix(data1,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_out = [data1; data1(end,:,:)];
 elseif dim == 4
    data1 = ncread(filedir, variable);

    latitude = [-90,90];
    longitude = [-180,180];

    data1 = rearrange_matrix(data1,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_out = [data1; data1(end,:,:,:)];
 else % dim = 2 
    data1 = ncread(filedir, variable);

    latitude = [-90,90];
    longitude = [-180,180];

    data1 = rearrange_matrix(data1,row,dim);

    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_out = [data1; data1(end,:)];
 end
    end