clear all
close all
% Plot the evolution within a cell over time

% Case List:
% - casenowaves := waves in ice module turned off
% - casezero := WIM on, with IFD being reset to 0 when ice conc. < puny
% - casenozero := WIM on, commenting out the reset

%% Preamble
 latit = -65; % desired latitude of cell -90 is south pole
 longi = 180; % desired longitude of cell 
 plot_title = 'Transect';
 cases = "8month";
 date = '2005-01-';
 times = ["01","02","03","05","06","07"];%["03600","07200","10800","14400","18000","21600","25200"];
 datapoints = length(times);
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'a1724548'; %a1724548, noahday, Noah
 variable = ["wave_sig_ht","peak_period","fsdrad","aice"];%,"fsdrad"]; %["aice","frazil"];%["wave_sig_ht","peak_period","fsdrad"];
 variable_label = ["Significant wave height (m)", "Peak period (s)", "Floe radius (m)", "Fractional ice concentration"];
 map_type = 'eqaazim'; %cassini
 no_cases = 1;
 temp = sprintf("Transect along %d degrees longitude", longi);
 plot_title = temp;
 filedir = [];

 no_variables = length(variable);
 t=tiledlayout(no_variables,1);
 
for k = 1:no_variables
 %% Data wrangling
 % a) Reading in netcdf
 filedir = [];
 for j = 1:length(times)
     if timestep == '1'
            filedir = [filedir; strcat('cases/', cases,'/history/iceh_inst.',date,'-',times(j),'.nc')];
     else
            filedir = [filedir; strcat('cases/',cases,'/history/iceh.',date,times(j),'.nc')];
     end

 % b) Find cell position
     [lat,lon,row] = grid_read(grid);
     [lat_out,lon_out] = lat_lon_finder(latit,longi,lat,lon);

 % c) Extracting variable data
     data(:,:) = data_format(filedir(j),variable(k),row,lat,lon);
     final_data(j) = data(lon_out,lat_out);
 end
 %% Plotting
 % Plot the variables for that cell
 x_axis = linspace(0,7,length(final_data));
 nexttile
 area(x_axis,final_data)
 
 xlabel('Time (days)')
 ylabel(variable_label(k), 'Interpreter','none')
 
 %nexttile
 %[len,wid] = size(data);
 %dummy_data = zeros(len,wid);
 %dummy_data(lon_out,lat_out) = data(lon_out,lat_out);
 %mapmakercell(dummy_data,lat,lon," ",map_type)
end
plot_title = sprintf("Cell at (%d S,%d E)",90-lat_out,lon_out);
plot_title = strcat(plot_title," (lat/lon)");
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
