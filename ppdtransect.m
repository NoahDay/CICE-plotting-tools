clear all
close all
%% Preamble
user = 'Noah'; %a1724548, noahday, Noah
case_name = 'swhdiag';
grid = 'gx1'; 
variable = 'peak_period'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir
variable_label = ["Peak period (s)"];
video_name = strcat('transect','_',variable, '_', case_name, '_', '2021_11_08', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 61;
day = 1;
month = "06";
year = "2005";
timestep = 'd';
no_cases = 1;
%% Plot Info
lon_transect = 1; % degrees East
 temp = sprintf("Transect along %d degrees longitude", lon_transect);
 plot_title = temp;
 filedir = "cases/ppdissue/history/iceh_01h.2005-06-03-03600.nc";
 %filedir = "cases/ppdissue/history/iceh_01h.2005-06-02-03600.nc";
 latit = 48; % 64 is equivalent to -90 to -30, 48 is approx -55
 latplot = 55;
 no_variables = length(variable);

%% Create dates and plot titles
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/20); % 0.5 = 2 seconds per frame
[lat,lon,row] = grid_read(grid);
 
%% Data wrangling
% a) create file directories

 
 
% b) Reading in netcdf and plot
  
   data(:,:) = data_format(filedir,'peak_period',row,lat,lon);
   [len,wid] = size(data);
    
   x_axis = linspace(latplot,79,latit);
   area(x_axis,data(lon_transect,latit:-1:1));

   xlabel('Latitude (degrees South)')
   ylabel(variable_label, 'Interpreter','none')
   limit = colorlims(variable);
   ylim(limit)
    


