close all
clear all
clc
addpath functions
% create video writer object
user = 'a1724548'; %a1724548, noahday, Noah
case_name = 'monthwim'; %ocnforcing
grid = 'gx1'; 
variable = "fsdrad"; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir, hi, uvel, vvel, Tair, frazil, iage
% afsd, dafsd_wave, Tsfc

video_name = strcat(variable, '_', case_name, '_', '2022_06_01', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
time_period = 'd'; %'1','d','m','y'
datapoints = 13;
%timestep = 7;
% Prestorm: 28/6
% Storm: 3/7
day = 6; month = 9; year = 2005;
sector = "SA";
if day < 9
    if month < 10
        date = sprintf('%d-0%d-0%d', year, month, day);
    else
        date = sprintf('%d-%d-%d', year, month, day);
    end
else
    if month < 10
        date = sprintf('%d-0%d-%d', year, month, day);
    else
        date = sprintf('%d-%d-%d', year, month, day);
    end
end
ssd = 1;
map_type = 'eqaazim'; %cassini
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/20); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);
ticker = 1;
SIC = 0.15; 

%Plotting

for i = 1:datapoints
    close all
   % Get the file name
    if ssd == 1
        filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',date,".nc");
    else
        filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
    end
   % Plot the map
   map_creator(filename, date, i, variable, grid, sector, user, SIC)
   % Update date
   %for j = 1:timestep
    date = update_date(date);
   %end
   if mod(i,floor(datapoints/10)) == 0
       clc
       fprintf('%g0%% complete\n',ticker);
       ticker = ticker + 1;
   end
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
 