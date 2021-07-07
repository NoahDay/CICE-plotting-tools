% create video writer object
 writerObj = VideoWriter('siconc.avi');
 plot_title = 'Sea Ice Concentration from 01 January 2005 - 10 January 2005';
 datapoints = 12;
 grid = 'gx3';
 user = 'noahday';
 variable = 'aice';
 % set the frame rate to one frame per second
 set(writerObj,'FrameRate',1);
 % open the writer
 open(writerObj);
 
 % sea ice area: 'aicen' or 'aice'
 % horizonal velocity: 'uvel'
 % vertical velocity: 'vvel'
 
if grid == 'gx3' % 3 degree grid
    for i = 1:datapoints
        filename = sprintf('case1y/history/iceh.2005-0%d.nc', i);
    if i > 9
        filename = sprintf('case1y/history/iceh.2005-%d.nc', i);
    end
    map_creator(filename, plot_title, i, variable, grid)
    end
elseif grid == 'gx1' % 1 degree grid
    for i = 1:datapoints
        filename = sprintf('casegx1/history/iceh.2005-01-0%d.nc', i);
    if i > 9
        filename = sprintf('casegx1/history/iceh.2005-01-%d.nc', i);
    end
    map_creator(filename, plot_title, i, 'aice', grid)
    end
end

    
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