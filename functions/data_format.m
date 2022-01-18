 function data_out = data_format(filedir,variable,row,lat,lon,dim)
 % Extract the data from a netCDF file.
 % filedir := file directory
 % variable := variable of data we want to extract
 % lat := grid latitude
 % lon := grid longitude
 
 % set number of dimensions to 2 by default
 if ~exist('dim', 'var')
    dim = 2; 
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