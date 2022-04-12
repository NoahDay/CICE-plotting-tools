% Search the data for a storm or cyclone.
% "The corrected ERA5 maximum wind speed was greater than 33 m/s" Vichi et
% al. (2019)
clear
close
addpath functions
addpath packages/quiverwcolorbar


user = 'a1724548'; %a1724548, noahday, Noah
case_name = 'ocnatmo';
grid = 'gx1'; 
time_period = 'd'; %'1','d','m','y'
datapoints = 365;
day = 1;
month = 1;
year = 2005;
sector = "SA";
if day > 9
    date = sprintf('%d-0%d-%d', year, month, day);
else
    date = sprintf('%d-0%d-0%d', year, month, day);
end

ticker = 1;
fig_count = 1;

%% Read in data
%ssd_dir = '/Volumes/Noah_SSD/run_data';
filedir = strcat('/Volumes/NoahDay5TB/','cases/',case_name);

for i = 1:datapoints
   % Get the file name
   filename = strcat(filedir,"/history/iceh.",date(i,:),".nc");
   % Find the maximum wind speed in that sector
   [wind_max(i) wind_mean(i) wind_min(i)] = storm_finder(filename, sector);
   
   % Update date
   date(i+1,:) = update_date(date(i,:));
    if mod(i,floor(datapoints/10)) == 0
       clc
       fprintf('%g0%% complete\n',ticker);
       ticker = ticker + 1;
    end
end    
        
storm_threshold = 24.6933;
cyclone_threshold =  32.9244;
idx_storm = wind_max > storm_threshold;
idx_cyclone = wind_max > cyclone_threshold;
% Winds greater than 24.6933 m/s are storms and greater than 32.9244 m/s
% are considered storms
% (http://www.bom.gov.au/marine/knowledge-centre/reference/wind.shtml).
    
storm_dates = date(idx_storm,:);
cyclone_dates = date(idx_cyclone,:);


%% Plot the ice and wind velocities every day that there is a storm
[lat,lon,row] = grid_read(grid);
[wid,len] = size(lat);

dim = 2;

[~, sector_mask] = data_format_sector(filename,"uatm",sector);
lat_sector = lat(sector_mask);
lon_sector = lon(sector_mask);
        
lat_vec = reshape(lat_sector,1,[]);
lon_vec = reshape(lon_sector,1,[]);

atm_var1 = "uatm"; % atm velocity (x)
atm_var2 = "vatm"; % atm velocity (y)
ice_var1 = "uvel"; % ice velocity (x)
ice_var2 = "vvel"; % ice velocity (y)

[len, wid] = size(storm_dates);
for k=1:len
    figure(fig_count)
    t = tiledlayout(1,2);
    title(t,strcat("Storm on ",convertCharsToStrings(storm_dates(k,:))))
    xlabel(t,'Longitude')
    ylabel(t,'Latitude')
        filename = strcat(filedir,"/history/iceh.",storm_dates(k,:),".nc");
        
        aice = data_format_sector(filename,"aice",sector);
        SIC = 0.15;
        [lat_ice_edge_new, lon_ice_edge_new, edge] = find_ice_edge(aice,SIC,sector,lat,lon);
        
        atm_data_x = data_format_sector(filename,atm_var1,sector);
        atm_data_y = data_format_sector(filename,atm_var2,sector);
        atm_data_x = atm_data_x(sector_mask);
        atm_data_y = atm_data_y(sector_mask);

        deltalat_atm_vec = reshape(atm_data_x,1,[]);
        deltalon_atm_vec = reshape(atm_data_y,1,[]);

        
        mag_atm = sqrt(deltalat_atm_vec.^2+deltalon_atm_vec.^2);
        scale = 5*1/max(mag_atm);
        
       nexttile
        quiverwcolorbar(lon_vec',-lat_vec',deltalon_atm_vec',-deltalat_atm_vec',scale,'bounds',[min(mag_atm) max(mag_atm)]);
            a = colorbar;
            a.Label.String = 'Velocity (m/s)';
            xlabel('Longitude (degrees East)');
            ylabel('Latitude (degrees South)');
            title(strcat("Wind velocities on ",convertCharsToStrings(storm_dates(k,:))))

        ice_data_x = data_format_sector(filename,ice_var1,sector);
        ice_data_y = data_format_sector(filename,ice_var2,sector);
        ice_data_x = ice_data_x(sector_mask);
        ice_data_y = ice_data_y(sector_mask);

        deltalat_ice_vec = reshape(ice_data_x,1,[]);
        deltalon_ice_vec = reshape(ice_data_y,1,[]);

        mag_ice = sqrt(deltalat_ice_vec.^2+deltalon_ice_vec.^2);
        scale = 0.6*1/max(mag_ice);
            
        nexttile
        if max(mag_ice) == 0    
            quiverwcolorbar(lon_vec',-lat_vec',deltalon_ice_vec',-deltalat_ice_vec',scale);
        else
            quiverwcolorbar(lon_vec',-lat_vec',deltalon_ice_vec',-deltalat_ice_vec',scale,'bounds',[min(mag_ice) max(mag_ice)]);
        end
            a = colorbar;
            a.Label.String = 'Velocity (m/s)';
            xlabel('Longitude (degrees East)');
            ylabel('Latitude (degrees South)');
            title(strcat("Ice velocities on ",convertCharsToStrings(storm_dates(k,:))))
        fig_count = fig_count + 1;    
end

%% Time series of wind speeds and frequencies of stormy days
end_date = date(end,:);
end_year = str2double(end_date(1:4));
end_month = str2double(end_date(6:7));
end_day = str2double(end_date(9:10));

% Storm time series
t1 = datetime(year,month,day);
t2 = datetime(end_year,end_month,end_day);
dates = datevec(t1:t2);

ts1 = timeseries(idx_storm,1:length(idx_storm));
ts1.Name = 'Storm frequency in the SA sector';
ts1.TimeInfo.Units = 'days';
month_name = datestr(datetime(1,month,1),'mmmm');
ts1.TimeInfo.StartDate = strcat(sprintf('%d',day),'-',month_name,'-',sprintf('%d',year));% Set start date.
ts1.TimeInfo.Format = 'dd mmm, yy';       % Set format for display on x-axis.

line_width = 3;

fig_count = fig_count + 1;
figure(fig_count)
t = tiledlayout(1,2);
nexttile
plot(ts1,'LineWidth',line_width)
ylim([0,2])

% Cyclone time series
t1 = datetime(year,month,day);
t2 = datetime(end_year,end_month,end_day);
dates = datevec(t1:t2);

ts2 = timeseries(idx_cyclone,1:length(idx_cyclone));
ts2.Name = 'Cyclone frequency in the SA sector';
ts2.TimeInfo.Units = 'days';
month_name = datestr(datetime(1,month,1),'mmmm');
ts2.TimeInfo.StartDate = strcat(sprintf('%d',day),'-',month_name,'-',sprintf('%d',year));     % Set start date.
ts2.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

line_width = 3;
nexttile
plot(ts2,'LineWidth',line_width)
ylim([0,2])


% Max wind speed
fig_count = fig_count + 1;
figure(fig_count)
t1 = datetime(year,month,day);
t2 = datetime(end_year,end_month,end_day);
dates = datevec(t1:t2);

ts2 = timeseries(wind_max,1:length(wind_max));
ts2.Name = 'Max wind speed (m/s)';
ts2.TimeInfo.Units = 'days';
month_name = datestr(datetime(1,month,1),'mmmm');
ts2.TimeInfo.StartDate = strcat(sprintf('%d',day-1),'-',month_name,'-',sprintf('%d',year));     % Set start date.
ts2.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

line_width = 3;
plot(ts2,'LineWidth',line_width)
grid on
yline(storm_threshold,'--')
yline(cyclone_threshold,'--')
xline(5,'--')
txt = 'Start of storm';
h = text(4.8,10,txt);
 set(h,'Rotation',90);
 ylim([0,35])

txt = 'Start of pre-storm';
h = text(0.2,10,txt);
 set(h,'Rotation',90); 


txt = 'Storm threshold (Beaufort Wind Scale)';
text(0,storm_threshold+1,txt)

txt = 'Cyclone threshold (Beaufort Wind Scale)';
text(0,cyclone_threshold+1,txt)


% Mean wind speed
fig_count = fig_count + 1;
figure(fig_count)
t1 = datetime(year,month,day);
t2 = datetime(end_year,end_month,end_day);
dates = datevec(t1:t2);

ts2 = timeseries(wind_mean,1:length(wind_max));
ts2.Name = 'Mean wind speed (m/s)';
ts2.TimeInfo.Units = 'days';
month_name = datestr(datetime(1,month,1),'mmmm');
ts2.TimeInfo.StartDate = strcat(sprintf('%d',day-1),'-',month_name,'-',sprintf('%d',year));     % Set start date.
ts2.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

line_width = 3;
plot(ts2,'LineWidth',line_width)
grid on
ylim([0,35])
yline(storm_threshold,'--')
yline(cyclone_threshold,'--')
xline(5,'--')
txt = 'Start of storm';
h = text(4.8,10,txt);
 set(h,'Rotation',90);

txt = 'Start of pre-storm';
h = text(0.2,10,txt);
 set(h,'Rotation',90); 


txt = 'Storm threshold (Beaufort Wind Scale)';
text(0,storm_threshold+1,txt)

txt = 'Cyclone threshold (Beaufort Wind Scale)';
text(0,cyclone_threshold+1,txt)


% Min wind speed
fig_count = fig_count + 1;
figure(fig_count)
t1 = datetime(year,month,day);
t2 = datetime(end_year,end_month,end_day);
dates = datevec(t1:t2);

ts2 = timeseries(wind_min,1:length(wind_max));
ts2.Name = 'Min wind speed (m/s)';
ts2.TimeInfo.Units = 'days';
month_name = datestr(datetime(1,month,1),'mmmm');
ts2.TimeInfo.StartDate = strcat(sprintf('%d',day-1),'-',month_name,'-',sprintf('%d',year));     % Set start date.
ts2.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

line_width = 3;
plot(ts2,'LineWidth',line_width)
 ylim([0,35])
grid on
yline(storm_threshold,'--')
yline(cyclone_threshold,'--')
xline(5,'--')
txt = 'Start of storm';
h = text(4.8,10,txt);
 set(h,'Rotation',90);

txt = 'Start of pre-storm';
h = text(0.2,10,txt);
 set(h,'Rotation',90); 


txt = 'Storm threshold (Beaufort Wind Scale)';
text(0,storm_threshold+1,txt)

txt = 'Cyclone threshold (Beaufort Wind Scale)';
text(0,cyclone_threshold+1,txt)


%% Investigate the stormiest day

% What day has the greatest winds
find(wind_max == max(wind_max));


% Calculate the ice edge
aice_data = data_format_sector(filename,"aice",sector);

lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);

SIC = 0.15;
[lat_ice_edge, lon_ice_edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);

 %% Functions
 