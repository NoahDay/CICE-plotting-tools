%% Ice drift
% Compare ice drift with wind speeds, ice concentration, internal stresses
% Also compare wind direction with ice direction.
clear
close all
addpath functions
addpath packages/bedmap2_toolbox_v4

location = [-65, 360-132.2];
SIC_threshold = 0.15;
%% Preamble
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'fixedwaves';

grid = 'gx1'; 
time_period = 'd'; %'1','d','m','y'
season = "winter";

day_init = 1;
month_init = 6;
year_init = 2006;
date = sprintf('%d-0%d-0%d', year_init, month_init, day_init);

dim = 2;
[lat,lon,row] = grid_read(grid);
filedir = strcat('cases/',case_name);
%% Read in Data
datapoints = 120;

[lat_out_north,lon_out] = lat_lon_finder(location(1),location(2),lat,lon); 
[lat_out_south,~] = lat_lon_finder(-80,location(2),lat,lon); 
for j = 1:datapoints
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    % Find the edge of the MIZ
    variable = "aice";
    aice_data = data_format(filename,variable,row,lat,lon,dim);
%     idx = 0;
%     for i = lat_out_north:-1:lat_out_south
%         idx = aice_data(lon_out,i) > SIC_threshold;
%         if idx == 1 % Found the ice edge
%             lat_out = i;
%         end
%     end
%     if i == lat_out_south
%         idx = 1;
%         lat_out = 1;
%     end
    lat_out = lat_out_north;
    % Wind speed
    variable = "uatm";
    atm_data_x = data_format(filename,variable,row,lat,lon,dim);
    variable = "vatm";
    atm_data_y = data_format(filename,variable,row,lat,lon,dim);

    atm_x = atm_data_x(lon_out,lat_out);
    atm_y = atm_data_y(lon_out,lat_out);

    wind_speed(j) = norm([atm_x,atm_y]);

    wind_direction(j) = atan(atm_y/atm_x);
    
    % Drift speed
    variable = "uvel";
    ice_data_x = data_format(filename,variable,row,lat,lon,dim);
    variable = "vvel";
    ice_data_y = data_format(filename,variable,row,lat,lon,dim);

    ice_x = ice_data_x(lon_out,lat_out);
    ice_y = ice_data_y(lon_out,lat_out);

    ice_drift_speed(j) = norm([ice_x,ice_y]);

    ice_direction(j) = atan(ice_y/ice_x);

    % Ice concentration
    aice(j) = aice_data(lon_out,lat_out);
% 
    % Internal stresses
    variable = "sigP";
    sigP_data = data_format(filename,variable,row,lat,lon,dim);
    sigP(j) = sigP_data(lon_out,lat_out);
    
    % FSD radius
    variable = "fsdrad";
    fsdrad_data = data_format(filename,variable,row,lat,lon,dim);
    fsdrad(j) = fsdrad_data(lon_out,lat_out);

    variable = "wave_sig_ht";
    swh_data = data_format(filename,variable,row,lat,lon,dim);
    swh(j) = swh_data(lon_out,lat_out);


    date = update_date(date);
end
%% Plotting
% Wind speed
t1 = datetime(year_init,month_init,day_init);
t2 = datetime(str2num(date(1:4)),str2num(date(6:7)),str2num(date(9:10)));
dates = datevec(t1:t2);

ts1 = timeseries(wind_speed,1:datapoints);
ts1.Name = 'Forcing wind speed';
ts1.TimeInfo.Units = 'days';
ts1.TimeInfo.StartDate = char(t1);     % Set start date.
ts1.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.


ts2 = timeseries(ice_drift_speed,1:datapoints);
ts2.Name = 'Ice drift';
ts2.TimeInfo.Units = 'days';
ts2.TimeInfo.StartDate = char(t1);     % Set start date.
ts2.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.
% 
% 
ts3 = timeseries(aice,1:datapoints);
ts3.Name = 'Ice concentration';
ts3.TimeInfo.Units = 'days';
ts3.TimeInfo.StartDate = char(t1);     % Set start date.
ts3.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

ts4 = timeseries(sigP,1:datapoints);
ts4.Name = 'Internal stress';
ts4.TimeInfo.Units = 'days';
ts4.TimeInfo.StartDate = char(t1);     % Set start date.
ts4.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.
% 
 ts5 = timeseries(fsdrad,1:datapoints);
ts5.Name = 'FSD radius';
ts5.TimeInfo.Units = 'days';
ts5.TimeInfo.StartDate = char(t1);     % Set start date.
ts5.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

ts6 = timeseries(swh,1:datapoints);
ts6.Name = 'SWH';
ts6.TimeInfo.Units = 'days';
ts6.TimeInfo.StartDate = char(t1);     % Set start date.
ts6.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

 ts5 = timeseries(wind_direction,1:datapoints);
ts5.Name = 'Wind direction';
ts5.TimeInfo.Units = 'days';
ts5.TimeInfo.StartDate = char(t1);     % Set start date.
ts5.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

ts6 = timeseries(ice_direction,1:datapoints);
ts6.Name = 'Ice drift direction';
ts6.TimeInfo.Units = 'days';
ts6.TimeInfo.StartDate = char(t1);     % Set start date.
ts6.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

ts7 = timeseries(wind_speed - ice_drift_speed,1:datapoints);
ts7.Name = 'Wind speed - Ice speed';
ts7.TimeInfo.Units = 'days';
ts7.TimeInfo.StartDate = char(t1);     % Set start date.
ts7.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.

ts8 = timeseries(wind_direction - ice_direction,1:datapoints);
ts8.Name = 'Wind direction - Ice direction';
ts8.TimeInfo.Units = 'days';
ts8.TimeInfo.StartDate = char(t1);     % Set start date.
ts8.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.




%% Plotting
line_width = 2;
t = tiledlayout(1,1);
 nexttile
% % Wind speed
%  plot(ts1,'LineWidth',line_width)
%  nexttile
% % Ice drift
%  plot(ts2,'LineWidth',line_width)
% nexttile
% % Ice concentration
% plot(ts3,'LineWidth',line_width)
%  nexttile
% % Internal stress
%  plot(ts4,'LineWidth',line_width)
% %nexttile
% % FSD radius
% %plot(ts5,'LineWidth',line_width)
% %nexttile
% % SWH
% %plot(ts6,'LineWidth',line_width)

%plot(ts7,'LineWidth',line_width)
%nexttile
plot(ts8,'LineWidth',line_width)



