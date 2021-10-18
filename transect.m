clear all
close all
% Plotting a transect of the CICE output

% Case List:
% - casenowaves := waves in ice module turned off
% - casezero := WIM on, with IFD being reset to 0 when ice conc. < puny
% - casenozero := WIM on, commenting out the reset

% Comparing CICE results with and without waves
%% Preamble
 plot_title = 'Transect';
 cases = "casecawcr";
 date = '2005-01-07';
 datapoints = 1;
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'noahday'; %a1724548, noahday, Noah
 variable = ["wave_sig_ht_d","peak_period_d","fsdrad_d","aice_d"];%,"fsdrad"]; %["aice","frazil"];%["wave_sig_ht","peak_period","fsdrad"];
 variable_label = ["Significant wave height (m)", "Peak period (s)", "Floe radius (m)", "Fractional ice concentration"];
 map_type = 'eqaazim'; %cassini
 no_cases = 1;
 lon_transect = 180;
 temp = sprintf("Transect along %d degrees longitude", lon_transect);
 plot_title = temp;
 filedir = [];
 latit = 64; % 64 is equivalent to -90 to -60
 no_variables = length(variable);

 %% Data wrangling
 % a) Reading in netcdf
 if timestep == '1'
    for i = 1:no_cases
        filedir = [filedir; strcat('cases/', cases(i),'/history/iceh_inst.',date,'.nc')];
    end
else
    for i = 1:no_cases
        filedir = [filedir; strcat('cases/',cases(i),'/history/iceh.',date,'.nc')];
    end
 end
% b) Extracting variable data
t=tiledlayout(no_variables,1);
for i = 1:no_variables
 [lat,lon,row] = grid_read(grid);
 data(:,:,i) = data_format(filedir,variable(i),row,lat,lon);
 [len,wid] = size(data);
 %dummy_data = zeros(len,wid);
 %dummy_data(37,1:latit) = data(37,1:latit,i);
 
 %% Plotting
 %mapmaker(dummy_data,lat,lon,plot_title,map_type)
 x_axis = linspace(-60,-90,latit);
 nexttile
 area(x_axis,data(lon_transect,latit:-1:1,i))
 
 xlabel('Latitude (degrees)')
 ylabel(variable_label(i), 'Interpreter','none')
 %text(-90,0,'South pole')
end
 title(t,plot_title, 'Interpreter','none')
