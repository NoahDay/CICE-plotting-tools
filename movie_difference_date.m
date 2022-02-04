close all
clear all
addpath functions
% create video writer object
user = 'noahday'; %a1724548, noahday, Noah
case_name = '2005_2006';
grid = 'gx1'; 
variable = 'aice'; 
video_name = strcat('difference_',variable, '_', case_name, '_', '2022_02_02_world', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 30;%242;
manual_limits = [-500,500];
% Date 1
day = 1;
month = 9;
year = 2005;
date(1,:) = sprintf('%d-0%d-0%d', year, month, day);
% Date 2
day = 1;
month = 9;
year = 2006;
date(2,:) = sprintf('%d-0%d-0%d', year, month, day);

sector = "world";

map_type = 'eqaazim'; %cassini
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/10); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);
ticker = 1;
plot_title = 'test'; 
%% Plotting

case_name = strcat('cases/',case_name);

for i = 1:datapoints
   % Get the file name
   filename1 = strcat(case_name,"/history/iceh.", date(1,:),".nc");
   filename2 = strcat(case_name,"/history/iceh.", date(2,:),".nc");
   % Plot the map
   plot_title = strcat("Difference between ", convertCharsToStrings(date(1,:)), " and ", convertCharsToStrings(date(2,:)));
   difference_creator(filename1, filename2, plot_title, i, variable, grid, sector, user, manual_limits)
   % Update date
   date(1,:) = update_date(date(1,:));
   date(2,:) = update_date(date(2,:));
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
 