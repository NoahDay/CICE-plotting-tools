clear all
close all
% Plot the aggreate floe size histogram

% Case List:
% - casenowaves := waves in ice module turned off
% - casezero := WIM on, with IFD being reset to 0 when ice conc. < puny
% - casenozero := WIM on, commenting out the reset

%% Preamble 
 plot_title = 'Histogram';
 cases = "8month";
 date = '2005-01-';
 times = ["01","02","03","05","06","07"];%["03600","07200","10800","14400","18000","21600","25200"];
 datapoints = length(times);
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'a1724548'; %a1724548, noahday, Noah
 variable = ["afsd"];
 variable_label = ["Areal floe size distribution"];
 map_type = 'eqaazim'; %cassini
 no_cases = 1;
 plot_title = "Aggregate plot";
 filedir = [];

 no_variables = length(variable);
 t=tiledlayout(datapoints,1);
 
for k = 1:datapoints
 %% Data wrangling
 % a) Reading in netcdf

    filedir = [filedir; strcat('cases/',cases,'/history/iceh.',date,times(k),'.nc')];


 % b) Find cell position
     [lat,lon,row] = grid_read(grid);
%     [lat_out,lon_out] = lat_lon_finder(latit,longi,lat,lon);

 % c) Extracting variable data
     data(:,:,:) = data_format(filedir(1),variable(1),row,lat,lon,3);
     %final_data(1,:) = mean(data(:,:,:),3);
     [len, wid, dep] = size(data); 
     summed = zeros(1,dep);
     count = 0;
     for i = 1:len
         for j = 1:wid
             if  ~all(isnan(data(i,j,:)))
                 for l = 1:dep
                    summed(1,l) = summed(1,l) + data(i,j,l);    
                 end
                 count = count + 1;
             end
         end
     end
     final_data = summed/count;
 %end
 %% Plotting
 % Plot the variables for that cell
 x_axis = [1:dep];% linspace(0,dep,length(final_data));
 nexttile
 hist(final_data,x_axis)
 
 if k == floor(datapoints/2)
     label_h = ylabel('Density', 'Interpreter','none', 'FontSize', 12);
 label_h.Position(2) =0; % change vertical position of ylabel
 end
 
end
 xlabel('Floe radius (left side of intervals)')
 
 %nexttile
 %[len,wid] = size(data);
 %dummy_data = zeros(len,wid);
 %dummy_data(lon_out,lat_out) = data(lon_out,lat_out);
 %mapmakercell(dummy_data,lat,lon," ",map_type)

%plot_title = sprintf("Cell at (%d S,%d E)",90-lat_out,lon_out);
plot_title = strcat("Aggregate of floe size radius globally");
title(t,plot_title, 'Interpreter','none')

 function map_out = mapmakercell(data,lat,lon,plot_title,map_type)
        color_map = seaicecolormap();
        if map_type == 'cassini'
            x_origin = -20;
        else
            x_origin = -90;
        end
        w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [x_origin 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 30);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        title(plot_title, 'Units', 'normalized', 'Position', [-0.0, 0.5, 0], 'Rotation',90,'Interpreter', 'none')%title(plot_title, 'Interpreter', 'none',[0.5, -0.1, 0])
        colorbar
    end
