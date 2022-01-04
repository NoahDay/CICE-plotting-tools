clear all
close all
% Plotting the FSD hsitograms along a specified transect of the CICE output
% Comparing CICE results with and without waves
%% Preamble
plot_title = 'FSD histograms at';
cases = ["swhdiag", "2monthbase"];
date = '2005-06-01';
datapoints = 60;
grid = 'gx1';
timestep = 'd'; % '1', 'd', 'm', 'y'
user = 'noahday'; %a1724548, noahday, Noah
no_cases = 1;
lon_transect = 130;
filedir = [];
latit = 48; % 64 is equivalent to -90 to -30, 48 is approx -55
lati = -65; % desired latitude
longi = 130; % desired longitide
latplot = 55;
no_variables = 2;%length(variable);
font_size = 14;
line_thick = 1;
time_period = 'd';
day = 1;
month = 6;
year = 2005;
 %% Setup movie
video_name = strcat('FSD', '_', 'cell_evolution', '_', '2021_11_15', '.avi');
writerObj = VideoWriter(video_name);
% set the frame rate to one frame per second
set(writerObj,'FrameRate',datapoints/10); % 0.5 = 2 seconds per frame
% open the writer
open(writerObj);
ticker = 1;
%% Find the longitude and latitude indeces
[lat,lon,row] = grid_read(grid);
[lat_out,lon_transect] = lat_lon_finder(lati,longi,lat,lon);
%% Plot titles
if time_period == 'd'
    plot_date = [];
    for p = 1:datapoints
        if ticker>30
            month = month + 1;
            ticker = 1;
        end
        day_temp = day-1+ticker;
        if day_temp < 10
            date_temp = sprintf("%d-%d-0%d", year, month, day_temp);
        else
            date_temp = sprintf('%d-%d-%d', year, month, day_temp);
        end
        date = datestr(datetime(date_temp));
        plot_date = [plot_date; date];
        ticker = ticker + 1;
    end
elseif time_period == 'm'
   plot_date = ["January 2005", "February  2005", "March  2005", "April  2005", "May 2005", "June 2005", "July 2005", "August 2005", "September 2005", "October 2005", "November 2005", "December 2005"]; 
elseif time_period == 'y'
   plot_date = ["January 2005", "January 2006", "January 2007", "January 2008", "January 2009", "January 2010"];
else
    plot_date = ["Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005","Jan 1 2005"];
end

%% File directories and plotting frames
    
for i = 1:datapoints
        file = [];
    for j =1:2
        ticker = 1;
        case_name = strcat('cases/',cases(j));
        if i < 31
            if i < 10    
               if time_period == 'd'
                   dir_name = strcat(case_name,'/history/iceh.2005-06-0%d.nc');
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
                   dir_name = strcat(case_name,'/history/iceh.2005-06-%d.nc');
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
        end
        
        if i > 30 % next month
            if ticker < 10    
               if time_period == 'd'
                   dir_name = strcat(case_name,'/history/iceh.2005-07-0%d.nc');
                   filename = sprintf(dir_name, ticker);
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
            if ticker > 9
               if time_period == 'd'
                   dir_name = strcat(case_name,'/history/iceh.2005-07-%d.nc');
                   filename = sprintf(dir_name, ticker);
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
            ticker = ticker + 1;
        end
        file = [file, filename];
    end
    t=tiledlayout(length(cases),1);
    % Now plot frames
    for k = 1:length(cases)
        data = data_format(file(k),'afsd',row,lat,lon,3);
        long_data_hist(1,:) =  data(lon_transect,lat_out,:);
        idx = isnan(long_data_hist);
        long_data_hist(idx) = 0;
        lat_lab = lat(lon_transect,lat_out);
        nexttile
        histogram('BinEdges',[0:12],'BinCounts',long_data_hist);
        if k == 1
                ylabel('With waves')
            else
                ylabel('Without waves')
        end
    end
    date = plot_date(i,:);
    temp_title = sprintf(' (%d,%d) (lat,lon), ',lati,longi);
    title(t,strcat(plot_title, temp_title,' ',date))
    if j == 2
        xlabel('Floe size categories')
    end
    figname = sprintf('image%d.png', i); 
    text = sprintf('/Users/%s/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames', user);
    fname = strcat('/Volumes/Noah SSD/MATLAB/PhD Project/CICE Plotting/frames','/',cases(1),'_',cases(2)); 
    saveas(gcf,fullfile(fname, figname));
end

%% Make the movie
 for k=1:datapoints
     % use imread to read the image
     fignamedir = strcat('frames/',cases(1),'_',cases(2),'/');
     filename = sprintf('image%d.png',k);
     figname = strcat(fignamedir,filename);
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