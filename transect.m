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
 cases = "8month";
 date = '2005-07-29';
 datapoints = 1;
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'noahday'; %a1724548, noahday, Noah
 ice_variable = ["fsdrad","aice","hi"];%,"fsdrad"]; %["aice","frazil"];%["wave_sig_ht","peak_period","fsdrad"];
 wave_variable = ["wave_sig_ht","peak_period", "ice thickness"];
 fsd_variable = ["dafsd_wave","dafsd_newi", "dafsd_weld"];
 variable_label_ice = ["Floe radius (m)", "Fractional ice concentration", "Ice thickness (m)"];
 variable_label_wave = ["Significant wave height (m)", "Peak period (s)"];
 variable_label_fsd = ["Change in FSD: waves","Change in FSD: new ice", "Change in FSD: welding"];
 map_type = 'eqaazim'; %cassini
 no_cases = 1;
 lon_transect = 130;
 temp = sprintf("Transect along %d E, ", lon_transect);
 temp = strcat(temp, date);
 plot_title = temp;
 filedir = [];
 latit = 48; % 64 is equivalent to -90 to -30, 48 is approx -55
 latplot = 55;
 no_variables = 2;%length(variable);
 font_size = 14;
 line_thick = 1;

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
 x_axis = linspace(latplot,79,latit);
%% Find the ice edge and land edge
[lat,lon,row] = grid_read(grid);
ice_edge(:,:) = data_format(filedir,'aice',row,lat,lon);
long_ice_edge = ice_edge(lon_transect,latit:-1:1);
pos = find(long_ice_edge > 0.05);
ice_pos = x_axis(pos(1));

land_edge(:,:) = data_format(filedir,'tmask',row,lat,lon);
long_land_edge = land_edge(lon_transect,latit:-1:1);
pos = find(long_land_edge);
land_pos = x_axis(pos(end));

lat_edges = 0;

%% Plot Ice statistics: Hice, Aice, FSD
%t=tiledlayout(2,1);
%nexttile
for i = 1:length(ice_variable)
 %[lat,lon,row] = grid_read(grid);
 data(:,:,i) = data_format(filedir,ice_variable(i),row,lat,lon);
 [len,wid] = size(data);
 long_data_ice(i,:) = data(lon_transect,latit:-1:1,i);
 idx = isnan(long_data_ice(i,:));
 long_data_ice(i,idx) = 0;
end
x_axis = linspace(latplot,79,latit);

figure(1)
% FSD Rad
ax1 = axes;
yyaxis left
plot(x_axis,long_data_ice(1,:),'LineWidth',line_thick)
pause(0.1)
xlabel('Latitude (degrees South)')
ylabel(variable_label_ice(1), 'Interpreter','none')
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
%grid(ax1,'on')
ytick = ax1.YTick;  
set(gca,'FontSize',font_size)
% aice
yyaxis right             
plot(x_axis,long_data_ice(2,:),'LineWidth',line_thick)
set(gca,'FontSize',font_size)

 
% xlabel('Latitude (degrees South)')
 ylabel(variable_label_ice(2), 'Interpreter','none')
% limit = colorlims(ice_variable(2));
% ylim([0,1.1*max(long_data_ice)]);

 %text(-90,0,'South pole')
 
ax2 = axes('position', ax1.Position);
plot(ax2,x_axis,long_data_ice(3,:), 'k','LineWidth',line_thick)
ylabel(variable_label_ice(3), 'Interpreter','none')
ytickformat('%.1f')
set(gca,'FontSize',font_size)
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
%grid(ax2,'on')
% Horizontally scale the y axis to alight the grid (again, be careful!)
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'             '}); 
 
 title(plot_title, 'Interpreter','none')
 xline(land_pos,'--',{'Land Edge'});
 xline(ice_pos,'--',{'Ice Edge'});

%% Plot the wave statistics: Hs, Tp
figure(2)
for i = 1:2
 %[lat,lon,row] = grid_read(grid);
 data(:,:,i) = data_format(filedir,wave_variable(i),row,lat,lon);
 [len,wid] = size(data);
 long_data_wave = data(lon_transect,latit:-1:1,i);
 idx = isnan(long_data_wave);
 long_data_wave(idx) = 0;
 %nexttile
 if rem(i, 2) == 0 % is even
    yyaxis right
 else
     yyaxis left
 end
 x_axis = linspace(latplot,79,latit);
 plot(x_axis,long_data_wave,'LineWidth',line_thick)
 %histogram('BinEdges',[54, x_axis], 'BinCounts', long_data_wave)
 xlabel('Latitude (degrees South)')
 ylabel(variable_label_wave(i), 'Interpreter','none')
 %limit = colorlims(wave_variable(i));
 ylim([0,1.1*max(long_data_wave)]);
 xline(land_pos,'--',{'Land Edge'});
 xline(ice_pos,'--',{'Ice Edge'});
 set(gca,'FontSize',font_size)
 %text(-90,0,'South pole')
end

%t=tiledlayout(2,1);
%nexttile
%% FSD changes
figure(3)
dim =3;
%t=tiledlayout(2,1);
%nexttile
for i = 1:length(fsd_variable)
 %[lat,lon,row] = grid_read(grid);
 data_read(:,:,:,i) = data_format(filedir,fsd_variable(i),row,lat,lon,dim);
 data(:,:,i) = data_read(:,:,1,i); % first category
 [len,wid] = size(data);
 long_data_fsd(i,:,:) = data(lon_transect,latit:-1:1,i);
 idx = isnan(long_data_fsd);
 long_data_fsd(idx) = 0;
 idx = isnan(long_data_ice(i,:));
 long_data_ice(i,idx) = 0;
end
x_axis = linspace(latplot,79,latit);

% FSD Rad
ax1 = axes;
yyaxis left
plot(x_axis,long_data_fsd(1,:),'LineWidth',line_thick)
pause(0.1)
xlabel('Latitude (degrees South)')
ylabel(variable_label_fsd(1), 'Interpreter','none')
ax1.XTickMode = 'manual'; 
ax1.YTickMode = 'manual'; 
ax1.YLim = [min(ax1.YTick), max(ax1.YTick)];  % see [4]
ax1.XLimMode = 'manual'; 
%grid(ax1,'on')
ytick = ax1.YTick;  
set(gca,'FontSize',font_size)
% aice
yyaxis right             
plot(x_axis,long_data_fsd(2,:),'LineWidth',line_thick)
set(gca,'FontSize',font_size)

 
% xlabel('Latitude (degrees South)')
 ylabel(variable_label_fsd(2), 'Interpreter','none')
% limit = colorlims(ice_variable(2));
% ylim([0,1.1*max(long_data_ice)]);

 %text(-90,0,'South pole')
 
ax2 = axes('position', ax1.Position);
plot(ax2,x_axis,long_data_fsd(3,:), 'k','LineWidth',line_thick)
ylabel(variable_label_fsd(3), 'Interpreter','none')
set(gca,'FontSize',font_size)
pause(0.1)                 % see [3]
ax2.Color = 'none'; 
%grid(ax2,'on')
% Horizontally scale the y axis to alight the grid (again, be careful!)
ax2.XLim = ax1.XLim; 
ax2.XTick = ax1.XTick; 
ax2.YLimMode = 'manual'; 
yl = ax2.YLim; 
ax2.YTick = linspace(yl(1), yl(2), length(ytick));      % see [2]
% horzontally offset y tick labels
ax2.YTickLabel = strcat(ax2.YTickLabel, {'         '}); 
 
 title(plot_title, 'Interpreter','none')
 xline(land_pos,'--',{'Land Edge'});
 xline(ice_pos,'--',{'Ice Edge'});
title(plot_title, 'Interpreter','none')
