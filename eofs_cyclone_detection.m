% Setup
clear all
clc
close all
addpath functions
addpath /Users/noahday/GitHub/CICE-analyser/processing

historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/2019/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end
%% Air pressure
filename.airPressure = '/Users/noahday/GitHub/CICE-plotting-tools/observations/jra55/psl_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201901010000-201912312100.nc';
airPressure = ncread(filename.airPressure, "psl");
lat_jra = ncread(filename.airPressure, "lat");
lon_jra = ncread(filename.airPressure, "lon");
time = ncread(filename.airPressure, "time");
sector = "AU";
coords = sector_coords(sector);
% September start iceh_inst.2007-07-01-75600.nc
jra_idx = 120*8; % May 1st
%181*8; % July 1st
[len,wid,~] = size(airPressure);
%%
model_timestep = 'd';

if model_timestep == 'h'
    for i = 1:ceil(n_files/3) % Data is every 3 hr
        for j = 1:len
            for k = 1:wid
                if lat_jra(k) > coords(2,1) && lat_jra(k) < coords(1,1) && mod(lon_jra(j)+180,360) > mod(coords(1,2)+180,360)&& mod(lon_jra(k)+180,360) < mod(coords(3,2)+180,360)
                    temp = squeeze(airPressure(j,k,jra_idx:jra_idx+ceil(n_files/3)-1));
                    airPressureInterp(j,k,1:n_files) = interp1(1:ceil(n_files/3),temp,linspace(1,ceil(n_files/3),n_files));
                else
                    airPressureInterp(j,k,1:n_files) = NaN;
                end
            end
        end
    end
elseif model_timestep == 'd'
    for i = 1:n_files % Data is every 8 times a day
        for j = 1:len
            for k = 1:wid
                if lat_jra(k) > coords(2,1) && lat_jra(k) < coords(1,1) && mod(lon_jra(j)+180,360) > mod(coords(1,2)+180,360)&& mod(lon_jra(k)+180,360) < mod(coords(3,2)+180,360)
                    temp = squeeze(airPressure(j,k,jra_idx+(i-1)*8));
                    temp_time = squeeze(airPressure(j,k,jra_idx+(i-1)*8));
                    %airPressureInterp(j,k,1:n_files) = interp1(1:ceil(n_files/3),temp,linspace(1,ceil(n_files/3),n_files));
                    airPressureInterp(j,k,i) = temp;
                else
                    airPressureInterp(j,k,i) = NaN;
                end
            end
        end
    end
end
clear temp
clear temp_time
%%
% Calculate mean air pressure
psl_mean = mean(airPressure,3);
imagescn(lon_jra,lat_jra,psl_mean)
cb = colorbar;

%% Calculate trend (hPa/yr)
psl_trend = 365.25*trend(airPressure,time,3);
imagescn(lon_jra,lat_jra,psl_trend)
cb = colorbar;
ylabel(cb,'Pressure trend')
cmocean('balance','pivot')

%% Eofs identify the modes of variability of a system

imagescn(lon_jra,lat_jra,eof(airPressure,1))
colorbar;
cmocean('balance','pivot')
title 'eof first mode'

%% Get time series at that location of interest
psl1 = squeeze(airPressure(2,50,:));
% (-62, 0.5)

plot(time,psl1)
datetick
% Deaseason the data
psl1_ds = deseason(psl1,time);

hold on
plot(time,psl1_ds)
hold off

%%
psl = airPressure(1:300,1:60,:)/1000;
[Lon,Lat] = meshgrid(lon_jra,lat_jra);
lon = lon_jra(1:300);
lat = lat_jra(1:50);
time = 1:2920;
psl_ds = deseason(psl,time);

imagescn(lon,lat,eof(psl_ds,1))
colorbar
cmocean('balance')



%%
psl_ds_dt = detrend3(psl_ds);

psl_anom_var = var(psl_ds_dt,[],3); % alont the third dimension

imagescn(lon,lat,psl_anom_var)
colorbar
%caxis([0,1])

%% What are the modes of variability in this deseasoned detrended data
lon = lon_jra;
lat = lat_jra;

% Calculate eofs
[eof_maps,pc,expv] = eof(psl_ds_dt,6);
% PC is the principal component time series
clf
subplot(3,2,1)
imagescn(lon,lat,eof_maps(:,:,1))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,2)
plot(time,pc(1,:))
axis tight
box off
datetick

subplot(3,2,3)
imagescn(lon,lat,eof_maps(:,:,2))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,4)
plot(time,pc(2,:))
axis tight
box off
datetick

subplot(3,2,5)
imagescn(lon,lat,eof_maps(:,:,3))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,6)
plot(time,pc(3,:))
axis tight
box off
datetick

sgtitle 'The first three principal components'


%%

h2.CData = eof_maps(:,:,1)*pc(1,1) + ...
           eof_maps(:,:,2)*pc(2,1) + ...
           eof_maps(:,:,3)*pc(3,1) + ...
           eof_maps(:,:,4)*pc(4,1) + ...
           eof_maps(:,:,5)*pc(5,1) + ...
           eof_maps(:,:,6)*pc(6,1);
% Reconstruct sst anomalies from first 5 modes
psl_ds_dt_r = reof(eof_maps,pc,1:6);
%t = datestr(time);
time_num = 1:2920;
for k = 1:120
    h1.CData = psl_ds_dt(:,:,k);
    h2.CData = psl_ds_dt_r(:,:,k);

    subplot(1,2,1)
    h1 = imagescn(lon,lat,h1.CData );
    title 'observed psl anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    
    subplot(1,2,2)
    h1 = imagescn(lon,lat,h2.CData);
    title 'reconstructed psl anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    pause(0.1)
    sgtitle(datestr(time(k),'dd/mm/yy'))
     if k == 1
        gif('reconstructed.gif','DelayTime',20/n_files,'resolution',100,'overwrite',true)
    else
        gif
    end
end



%%
clf
coord = [10,40];
plot(time_num,squeeze(psl_ds_dt(coord(1),coord(2),:)))
hold on
plot(time_num,squeeze(psl_ds_dt_r(coord(1),coord(2),:)))
hold off
legend({'Obs','PC1'})



%% PCA of CICE data
% aice
[lat,lon] = grid_read("om2");
time_cice = 1:length(filenames);
for i = 1:length(filenames)
    daidtd(:,:,i) = data_format_sector(filenames(i,:),"aice","SH");
    %daidtt(:,:,i) = data_format_sector(filenames(i,:),"daidtt","SH");
    daidtd_sector(:,:,i) = daidtd(:,1:50,i);
    %daidtt_sector(:,:,i) = daidtt(:,1:50,i);
end
lon_sector = lon(:,1:50);
lat_sector = lat(:,1:50);
%% Deaseason the data
daidtd_ds = deseason(daidtd,time_cice);
daidtt_ds = deseason(daidtt,time_cice);


% Detrend
daidtd_ds_dt = detrend3(daidtd_ds);
daidtt_ds_dt = detrend3(daidtt_ds);
%psl_anom_var = var(psl_ds_dt,[],3); % along the third dimension
%%
imagescn(lon',lat',eof(daidtd_ds,1)')
colorbar;
cmocean('balance','pivot')
title 'eof first mode'

%% Calculate eofs
%[eof_maps,pc,expv] = eof(daidtd_sector,6);
clf
for i = 1:6
    subplot(6,1,i)
    plot(pc(i,:))
end
%%
% Reconstruct daidtd anomalies from first 5 modes
daidtd_ds_dt_r = reof(eof_maps,pc,1:5);
%t = datestr(time);

for k = 1:120
    h1_CData = daidtd_sector(:,:,k);
    h2_CData = daidtd_ds_dt_r(:,:,k);

    subplot(1,2,1)
    h1 = imagescn(lon_sector',lat_sector',h1_CData');
    title 'observed daidtd anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    
    subplot(1,2,2)
    h2 = imagescn(lon_sector',lat_sector',h2_CData');
    title 'reconstructed daidtd anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    %pause(0.1)
    sgtitle(datestr(time_cice(k)))%,'dd/mm/yy'))
     if k == 1
        gif('reconstructed_daidtd.gif','DelayTime',20/n_files,'resolution',100,'overwrite',true)
    else
        gif
    end
end
%%
reconstructed_eof(filenames,"wave_sig_ht")

%%
pca_cube = load('pc_cube_june_sep_2019.mat');
data = pca_cube.pc_cube;
[~,~,dep] = size(data);
[lat,lon] = grid_read("om2");
    
time_cice = 1:dep;
lon_sector = lon(:,1:60);
lat_sector = lat(:,1:60);
% Deaseason the data
%daidtd_ds = deseason(daidtd,time_cice);
%daidtt_ds = deseason(daidtt,time_cice);

data_sector = data(:,1:60,:);
% Detrend
data_sector_dt = detrend3(data_sector);
data_sector_dt = detrend3(data_sector);
%psl_anom_var = var(psl_ds_dt,[],3); % along the third dimension
%
%imagescn(lon',lat',eof(daidtd_ds,1)')
%colorbar;
%cmocean('balance','pivot')
%title 'eof first mode'

% Calculate eofs
[eof_maps,pc,expv] = eof(data_sector_dt,6);

% Reconstruct daidtd anomalies from first 5 modes
data_r = reof(eof_maps,pc,1:5);
%t = datestr(time);
    
for k = 1:dep
    h1_CData = data_sector_dt(:,:,k);
    h2_CData = data_r(:,:,k);

    subplot(2,1,1)
    h1 = imagescn(lon_sector',lat_sector',h1_CData');
    title 'observed anomaly'
    colorbar
    cmocean bal
    caxis([-max(max(max(abs(data_sector_dt)))),max(max(max(abs(data_sector_dt))))])
    
    subplot(2,1,2)
    h2 = imagescn(lon_sector',lat_sector',h2_CData');
    title 'reconstructed anomaly'
    colorbar
    cmocean bal
    caxis([-max(max(max(abs(data_sector_dt)))),max(max(max(abs(data_sector_dt))))])
    %pause(0.1)
    DateVector = [2019,5,31,0,0,0];
    sgtitle(datestr(datetime(datestr(DateVector))+k))%,'dd/mm/yy'))
     if k == 1
        gif('pc1_anomalies.gif','DelayTime',0.3,'resolution',100,'overwrite',true)
    else
        gif
    end
end

%% FUNCTIONS

function reconstructed_eof(filenames,variable)
    [lat,lon] = grid_read("om2");
    time_cice = 1:length(filenames);
    for i = 1:length(filenames)
        if variable == "atmspd"
            udata(:,:,i) = data_format_sector(filenames(i,:),"uatm","SH");
            vdata(:,:,i) = data_format_sector(filenames(i,:),"vatm","SH");
            data(:,:,i) = sqrt(udata(:,:,i).^2 + vdata(:,:,i).^2);
        elseif variable == "vel"
            udata(:,:,i) = data_format_sector(filenames(i,:),"uvel","SH");
            vdata(:,:,i) = data_format_sector(filenames(i,:),"vvel","SH");
            data(:,:,i) = sqrt(udata(:,:,i).^2 + vdata(:,:,i).^2);
        else
            data(:,:,i) = data_format_sector(filenames(i,:),variable,"SH");
        end
        data_sector(:,:,i) = data(:,1:60,i);
    end
    lon_sector = lon(:,1:60);
    lat_sector = lat(:,1:60);
    % Deaseason the data
    %daidtd_ds = deseason(daidtd,time_cice);
    %daidtt_ds = deseason(daidtt,time_cice);
    
    
    % Detrend
    data_sector_dt = detrend3(data_sector);
    data_sector_dt = detrend3(data_sector);
    %psl_anom_var = var(psl_ds_dt,[],3); % along the third dimension
    %
    %imagescn(lon',lat',eof(daidtd_ds,1)')
    %colorbar;
    %cmocean('balance','pivot')
    %title 'eof first mode'
    
    % Calculate eofs
    [eof_maps,pc,expv] = eof(data_sector_dt,6);
    
    % Reconstruct daidtd anomalies from first 5 modes
    data_r = reof(eof_maps,pc,1:5);
    %t = datestr(time);
    
    for k = 1:120
        h1_CData = data_sector_dt(:,:,k);
        h2_CData = data_r(:,:,k);
    
        subplot(2,1,1)
        h1 = imagescn(lon_sector',lat_sector',h1_CData');
        title 'observed anomaly'
        cmocean bal
        caxis([-max(max(max(abs(data_sector_dt)))),max(max(max(abs(data_sector_dt))))])
        
        subplot(2,1,2)
        h2 = imagescn(lon_sector',lat_sector',h2_CData');
        title 'reconstructed anomaly'
        cmocean bal
       caxis([-max(max(max(abs(data_sector_dt)))),max(max(max(abs(data_sector_dt))))])
        %pause(0.1)
        DateVector = [2019,4,30,0,0,0];
        sgtitle(datestr(datetime(datestr(DateVector))+k))%,'dd/mm/yy'))
         if k == 1
            gif_name = strcat('reconstructed_',convertStringsToChars(variable),'.gif');
            gif(gif_name,'DelayTime',20/length(filenames),'resolution',100,'overwrite',true)
        else
            gif
        end
    end

end
























