close all
clear all
addpath functions
% create video writer object
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'fixedwaves';
grid = 'gx1'; 
variable = 'aice'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir, hi, uvel, vvel
video_name = strcat(variable, '_', case_name, '_', '2022_02_10', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 10;%242;
day = 1;
month = 2;
year = 2006;
sector = "world";
date = sprintf('%d-0%d-0%d', year, month, day);
map_type = 'eqaazim'; %cassini
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/10); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);
ticker = 1;
 
%% Plotting
user = 'noahday'; %a1724548, noahday, Noah
%ssd_dir = '/Users/noahday/Maths1/';%'/Volumes/Noah_SSD/run_data';

for i = 1:datapoints
   % Get the file name
   filename = strcat('cases/',case_name,"/history/iceh.",date,".nc");
   % Plot the map
   map_creator(filename, date, i, variable, grid, sector, user)
   % Update date
   date = update_date(date);
end    
        

    
% iterate over each image
 for k=1:datapoints
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     % resize the imagei_vec
     %img = imresize(img,2);
     % convert the image to a frame using im2frame
     frame = im2frame(img);
     % write the frame to the video
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);
 
 %% Functions
 