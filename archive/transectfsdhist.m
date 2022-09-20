clear all
close all
% Plotting the FSD hsitograms along a specified transect of the CICE output
% Comparing CICE results with and without waves
%% Preamble
 plot_title = 'Transect';
 cases = ["8month", "1yearnowaves"];
 date = '2005-07-29';
 datapoints = 1;
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'noahday'; %a1724548, noahday, Noah
 ice_variable = ["afsd","aice","hi"];%,"fsdrad"]; %["aice","frazil"];%["wave_sig_ht","peak_period","fsdrad"];
 wave_variable = ["wave_sig_ht","peak_period", "ice thickness"];
 variable_label_ice = ["Areal FSD", "Fractional ice concentration", "Ice thickness (m)"];
 variable_label_wave = ["Significant wave height (m)", "Peak period (s)"];
 map_type = 'eqaazim'; %cassini
 no_cases = 1;
 longi = 130;
 lati = -55;
 temp = sprintf("Transect along %d degrees longitude", longi);
 plot_title = temp;
 filedir = [];
 %latit = 48; % 64 is equivalent to -90 to -30, 48 is approx -55
 latplot = 55;
 no_variables = 2;%length(variable);
 font_size = 14;
 line_thick = 1;

 %% Plot CICE with and without Waves
 if timestep == '1'
    for i = 1:length(cases)
        filedir = [filedir; strcat('cases/', cases(i),'/history/iceh_inst.',date,'.nc')];
    end
else
    for i = 1:length(cases)
        filedir = [filedir; strcat('cases/',cases(i),'/history/iceh.',date,'.nc')];
    end
 end
 

%% Find the ice edge and land edge
[lat,lon,row] = grid_read(grid);
[latit,lon_transect] = lat_lon_finder(lati,longi,lat,lon);

x_axis = linspace(latplot,79,latit);
fsd_cat = ncread(filedir(1), 'NFSD');
fsd_range = [0; fsd_cat];
fsd_bins = [fsd_range(1:12),fsd_range(2:end)];

ice_edge(:,:) = data_format(filedir(1),'aice',row,lat,lon);
long_ice_edge = ice_edge(lon_transect,latit:-1:1);
pos = find(long_ice_edge > 0.05);
ice_pos = x_axis(pos(1));

land_edge(:,:) = data_format(filedir(1),'tmask',row,lat,lon);
long_land_edge = land_edge(lon_transect,latit:-1:1);
pos = find(long_land_edge);
land_pos = x_axis(pos(end));
transect_length = length(latit:-1:1);

%% Plot FSD 
t=tiledlayout(1,2)%(2,14);%transect_length);
for j=1:length(cases)
    data = data_format(filedir(j),'afsd',row,lat,lon,3);
    count = 1;
    for i = 1:transect_length
     if max(data(lon_transect,latit-i+1,:))>0
        long_data_hist(count,:) =  data(lon_transect,latit-i+1,:);
        lat_lab(count) = lat(lon_transect,latit-i+1);
        count = count + 1;
     end
    end
    [len,wid,dep] = size(long_data_hist);
    for i = 14%1:len
        nexttile
        histogram('BinEdges',[0:12],'BinCounts',long_data_hist(i,:));
        title(lat_lab(i))
        if i == 14
            if j == 1
                ylabel('With waves','Units', 'normalized', 'Position', [-0.15, 0.5, 0])
            else
                ylabel('Without waves','Units', 'normalized', 'Position', [-0.15, 0.5, 0])
            end
        end
        xlabel('FSD category')
    end
end
%xline(land_pos,'--',{'Land Edge'});
%xline(ice_pos,'--',{'Ice Edge'});
%title(t,sprintf('FSD histograms (12 categories) along %d degrees',longi))

