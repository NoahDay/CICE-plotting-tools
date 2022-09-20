clear all
close all
%% Preamble
user = 'Noah'; %a1724548, noahday, Noah
case_name = '8month';
grid = 'gx1'; 
variable = 'peak_period'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir
variable_label = ["Peak period (s)"];
video_name = strcat('transect','_',variable, '_', case_name, '_', '2021_11_08', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 242;
day = 1;
month = 1;
year = 2005;
date = sprintf('%d-0%d-0%d', year, month, day);
timestep = 'd';
no_cases = 1;
fontSize = 14; 
%% Plot Info
lon_transect = 1; % degrees East
 temp = sprintf("Transect along %d degrees longitude", lon_transect);
 plot_title = temp;
 filedir = [];
 latit = 48; % 64 is equivalent to -90 to -30, 48 is approx -55
 latplot = 55;
 no_variables = length(variable);

%% Create dates and plot titles
map_type = 'eqaazim'; %cassini
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/20); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);


[lat,lon,row] = grid_read(grid);
 
%% Plotting
% a) create file directories
case_name = strcat('cases/',case_name);
  
for i = 1:datapoints
    
   % Get the file name
   filename = strcat(case_name,"/history/iceh.",date,".nc");
   % Get data
   data = ncread(filename, variable);
   [~, ~, n] = size(data);
  
   data = data_format(filename,variable,row,lat,lon);
   [len,wid] = size(data);
    
   x_axis = linspace(latplot,79,latit);
   plot_data = data(lon_transect,latit:-1:1);
   idx = isnan(plot_data);
   plot_data(idx) = 0;
   area(x_axis,plot_data);

   xlabel('Latitude (degrees South)')
   ylabel(variable_label, 'Interpreter','none')
   limit = colorlims(variable);
   ylim(limit)
    

   title(date, 'FontSize', fontSize);
   figname = sprintf('image%d.png', i); 
   text = sprintf('/Users/%s/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames', user);
   fname = '/Volumes/Noah SSD/MATLAB/PhD Project/CICE Plotting/frames'; 
   saveas(gcf,fullfile(fname, figname));
    
    % Update date
    date = update_date(date);
 end
%% Make Movie
     % iterate over each image
 for k=1:datapoints
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     frame = im2frame(img);
     % write the frame to the video
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);