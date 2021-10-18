% create video writer object
user = 'a1724548'; %a1724548, noahday, Noah
case_name = 'withwavesnoerror';
grid = 'gx1'; 
variable = 'fsdrad'; % wave_sig_ht, peak_period, fsdrad, aice, mean_wave_dir
video_name = strcat(variable, '_', case_name, '_', '2021_10_18', '.avi');
writerObj = VideoWriter(video_name);
time_period = 'd'; %'1','d','m','y'
datapoints = 25;
 
 map_type = 'eqaazim'; %cassini
 % set the frame rate to one frame per second
 set(writerObj,'FrameRate',1); % 0.5 = 2 seconds per frame
 % open the writer
 open(writerObj);
 
 % sea ice area: 'aicen' or 'aice'
 % horizonal velocity: 'uvel'
 % vertical velocity: 'vvel'
 
if time_period == 'd'
   plot_title = ["01 January 2005", "02 January 2005", "03 January 2005", "04 January 2005", "05 January 2005", "06 January 2005", "07 January 2005", "08 January 2005", "09 January 2005", "10 January 2005", "11 January 2005", "12 January 2005", "13 January 2005", "14 January 2005", "15 January 2005", "16 January 2005", "17 January 2005", "18 January 2005","19 January 2005","20 January 2005","21 January 2005","22 January 2005","23 January 2005","24 January 2005","25 January 2005","26 January 2005","27 January 2005","28 January 2005","29 January 2005"];
elseif time_period == 'm'
   plot_title = ["January 2005", "February  2005", "March  2005", "April  2005", "May 2005", "June 2005", "July 2005", "August 2005", "September 2005", "October 2005", "November 2005", "December 2005"]; 
elseif time_period == 'y'
   plot_title = ["January 2005", "January 2006", "January 2007", "January 2008", "January 2009", "January 2010"];
else
    plot_title = ["Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005"];
end
 
case_name = strcat('cases/',case_name);
 
%if grid == 'gx3' % 3 degree grid
    for i = 1:datapoints
    if i < 10    
       if time_period == 'd'
           dir_name = strcat(case_name,'/history/iceh.2005-01-0%d.nc');
           filename = sprintf(dir_name, i);
       elseif time_period == 'm'
            dir_name = strcat(case_name,'/history/iceh.2005-0%d.nc');
            filename = sprintf(dir_name, i);
       elseif time_period == 'y'
            dir_name = strcat(case_name,'/history/iceh.200%d.nc');
            filename = sprintf(dir_name, i);
       else
           dir_name = strcat(case_name,'/history/iceh_inst.2005-01-01-0%d.nc');
            filename = sprintf(dir_name, i);
       end
    end
    if i > 9
       if time_period == 'd'
           dir_name = strcat(case_name,'/history/iceh.2005-01-%d.nc');
           filename = sprintf(dir_name, i);
       elseif time_period == 'm'
            dir_name = strcat(case_name,'/history/iceh.2005-%d.nc');
            filename = sprintf(dir_name, i);
       elseif time_period == 'y'
            dir_name = strcat(case_name,'/history/iceh.200%d.nc');
            filename = sprintf(dir_name, i);
       else
           dir_name = strcat(case_name,'/history/iceh_inst.2005-01-01-%d.nc');
            filename = sprintf(dir_name, i);
       end
    end
    map_creator(filename, plot_title, i, variable, grid, map_type, user)
    end
% elseif grid == 'gx1' % 1 degree grid
%     for i = 1:datapoints
%         filename = sprintf('caseppd/history/iceh.2005-0%d.nc', i);
%     if i > 9
%         filename = sprintf('caseppd/history/iceh.2005-%d.nc', i);
%     end
%     map_creator(filename, plot_title, i, variable, grid, map_type, user)
%     end
% end

    
     % iterate over each image
 for k=1:datapoints
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     % resize the imagei_vec
     %img = imresize(img,...);
     % convert the image to a frame using im2frame
     frame = im2frame(img);
     % write the frame to the video
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);