close all
clear all
% create video writer object
user = 'a1724548'; %a1724548, noahday, Noah
case_name = '1yearnowaves';%'8month';
grid = 'gx1'; 
variable = 'aice'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir, hi, uvel, vvel
video_name = strcat(variable, '_', case_name, '_', '2021_11_30', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 365;%242;
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
 
%% Plotting

case_name = strcat('cases/',case_name);

for i = 1:datapoints
   % Get the file name
   filename = strcat(case_name,"/history/iceh.",date,".nc");
   % Plot the map
   map_creator(filename, date, i, variable, grid, map_type, user)
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
 