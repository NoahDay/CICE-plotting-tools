clear all
close all
addpath functions
%% Preamble

close all
clear all
addpath functions
% create video writer object
user = 'noahday'; %a1724548, noahday, Noah
case_name = '8month';
grid = 'gx1'; 
variable = 'fsdrad'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir, hi, uvel, vvel
video_name = strcat(variable, '_', case_name, '_', '2021_11_30', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 180%365;%242;
day = 1;
month = 1;
year = 2005;
date = sprintf('%d-0%d-0%d', year, month, day);
map_type = 'eqaazim'; %cassini
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/10); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);
ticker = 1;

%% Loading data

% Grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');

dim = 2;
lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);

lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

color_map = seaicecolormap();

% Load ice shelf data
addpath packages/bedmap2_toolbox_v4

shelf = bedmap2_data('icemask');
shelf = 0.0*(shelf == 1); 
shelf(shelf==0) = NaN; % 1 iceshelf
[latshelf,lonshelf] = bedmap2_data('latlon');

% CICE data
ave_data = aggregate_data(case_name,date,datapoints,variable,dim);

%% Mapping

latitude = [-90,-30];
longitude = [-180,180];

w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,ave_data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    pcolorm(latshelf,lonshelf,shelf)  
    colorbar

    %title("Representative FSD per cell (m) on 31-08-2005",'interpreter','latex','FontSize', 18)

%% Functions

function ave_data = aggregate_data(case_name,date,datapoints,variable,dim)
    for i = 1:datapoints
        filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
        data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
        data1 = data;
        data = data(:,:,1);
        data = rearrange_matrix(data,37,dim);
        total_data(:,:,i) = [data; data(end,:)];
        
        % update date
        date = update_date(date);
    end
    ave_data = mean(total_data,3);
end