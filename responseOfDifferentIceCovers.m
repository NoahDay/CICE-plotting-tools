close all
clear all
% Add libraries
addpath functions
sector = "EA";
% Setup
historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/2019/may-sep/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end

% Read 
load('kmeans_may_sep_2019_4_classes.mat');
idx = kmeans_cluster.idx;
row_idx = kmeans_cluster.row_idx;
label_vec = kmeans_cluster.label_vec;
C = kmeans_cluster.C;
X_new = kmeans_cluster.X_new;
num_clusters = size(C,1);
%% Read in data cubes
for i = 1:n_files
    if i == 1
        NFSD = ncread(filenames(1,:),"NFSD");
        [floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
    end

    uatm(:,:,i) = data_format_sector(filenames(i,:),"uatm",sector);
    vatm(:,:,i) = data_format_sector(filenames(i,:),"vatm",sector);
    swh(:,:,i) = data_format_sector(filenames(i,:),"wave_sig_ht",sector);
    uvel(:,:,i) = data_format_sector(filenames(i,:),"uvel",sector);
    vvel(:,:,i) = data_format_sector(filenames(i,:),"vvel",sector);
    sst(:,:,i) = data_format_sector(filenames(i,:),"sst",sector);
    Tair(:,:,i) = data_format_sector(filenames(i,:),"Tair",sector);
    aice(:,:,i) = data_format_sector(filenames(i,:),"aice",sector);
    iage(:,:,i) = data_format_sector(filenames(i,:),"iage",sector);
    daidtd(:,:,i) = data_format_sector(filenames(i,:),"daidtd",sector);
    daidtt(:,:,i) = data_format_sector(filenames(i,:),"daidtt",sector);
    peak_period(:,:,i) = data_format_sector(filenames(i,:),"peak_period",sector);
    fsdrad(:,:,i) = data_format_sector(filenames(i,:),"fsdrad",sector);
    fsdrad(:,:,i) = data_format_sector(filenames(i,:),"fsdrad",sector);
    frzmlt(:,:,i) = data_format_sector(filenames(i,:),"frzmlt",sector); % 
    
    afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector); % m
    newi(:,:,:) = data_format_sector(filenames(i,:),"dafsd_newi",sector); % m
    weld(:,:,:) = data_format_sector(filenames(i,:),"dafsd_weld",sector); % m
    latg(:,:,:) = data_format_sector(filenames(i,:),"dafsd_latg",sector); % m
    latm(:,:,:) = data_format_sector(filenames(i,:),"dafsd_latm",sector); % m
    wave(:,:,:) = data_format_sector(filenames(i,:),"dafsd_wave",sector); % m
    for j = 1:numel(NFSD)
        itd1(:,:,j,i) = afsdn(:,:,j,1).*floe_binwidth(j);

        dafsd_newi(:,:,j,i) = newi(:,:,j).*floe_binwidth(j);
        dafsd_weld(:,:,j,i) = weld(:,:,j).*floe_binwidth(j);
        dafsd_latg(:,:,j,i) = latg(:,:,j).*floe_binwidth(j);
        dafsd_latm(:,:,j,i) = latm(:,:,j).*floe_binwidth(j);
        dafsd_wave(:,:,j,i) = wave(:,:,j).*floe_binwidth(j);
    end
end
row_file = [0,cumsum(row_idx)];
index_lat = X_new(:,end-1) == X_new(1,end-1);
index_lon = X_new(:,end) == X_new(1,end);
index_both = index_lat.*index_lon;
%%
[lat,lon] = grid_read('om2');
for file_number = 1:length(row_file)-1
    row_file = find(index_both);
    row_vec = row_file(file_number)+1:row_file(file_number+1);
    file_idx = idx(row_vec);
    [~, sector_mask] = data_format_sector(filenames(file_number,:),'aice',sector);
    X_map = [file_idx, X_new(row_vec,end-1:end)];
    [len,wid] = size(lat);
    k_means = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(lon,lat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos) = file_idx(i); 
        k_means(~sector_mask) = NaN;
    end
    k_means_cube(:,:,file_number) = k_means;
end

%% Get edge for the MIZ
clear swh_vec windspeed_vec winddir_vec icespeed_vec icedir_vec Tair_vec sst_vec daidtt_vec daidtd_vec aice_vec iage_vec peak_period_vec edge lat_vec
date_label = datetime(dirdates);
lon_pos = 10;
clust_vec = [2,4];
n_clust_vec = length(clust_vec);

longs = 80;
for lon_pos = longs
for n_clust = 1:n_clust_vec
    for i = 1:n_files-1
        [lat_ice_edge_new, lon_ice_edge_new, edge(i,:)] = find_edge(k_means_cube(:,:,i),clust_vec(n_clust),sector,lat,lon);
        %if n_clust == 2
        edge(i,:) = edge(i,:) - 1;
        %end
    end

    [~,wid] = size(edge);
    for i = 1:n_files-2
        % Difference in ice edge
        delta_edge(i,:) = edge(i+1,:) - edge(i,:);
        for j = 1:length(delta_edge(i,:))
            if delta_edge(i,j) > 0
                edge_idx(i,j) = 1;
            end
        end
    end
    
    data_pos = find(lon_pos == longs);
    for i = 1:length(edge(:,lon_pos)) % Each file
        if squeeze(aice(lon_pos,edge(i,lon_pos),i)) == 0
            %aice(lon_pos,edge(i,lon_pos),i) = NaN;
            disp('error')
        end
        lat_vec(i,data_pos,n_clust) = squeeze(lat(lon_pos,edge(i,lon_pos)));
        swh_vec(i,data_pos,n_clust) = squeeze(swh(lon_pos,edge(i,lon_pos),i));
        windspeed_vec(i,data_pos,n_clust) = sqrt(squeeze(uatm(lon_pos,edge(i,lon_pos),i)).^2 + squeeze(vatm(lon_pos,edge(i,lon_pos),i)).^2);
        winddir_vec(i,data_pos,n_clust) = atan(squeeze(vatm(lon_pos,edge(i,lon_pos),i))/squeeze(uatm(lon_pos,edge(i,lon_pos),i)));
        icespeed_vec(i,data_pos,n_clust) = sqrt(squeeze(uvel(lon_pos,edge(i,lon_pos),i)).^2 + squeeze(vvel(lon_pos,edge(i,lon_pos),i)).^2);
        icedir_vec(i,data_pos,n_clust) = atan(squeeze(vvel(lon_pos,edge(i,lon_pos),i))/squeeze(uvel(lon_pos,edge(i,lon_pos),i)));
        Tair_vec(i,data_pos,n_clust) = squeeze(Tair(lon_pos,edge(i,lon_pos),i));
        sst_vec(i,data_pos,n_clust) = squeeze(sst(lon_pos,edge(i,lon_pos),i));
        daidtt_vec(i,data_pos,n_clust) = squeeze(daidtt(lon_pos,edge(i,lon_pos),i));
        daidtd_vec(i,data_pos,n_clust) = squeeze(daidtd(lon_pos,edge(i,lon_pos),i));
        
        aice_vec(i,data_pos,n_clust) = squeeze(aice(lon_pos,edge(i,lon_pos),i));
        iage_vec(i,data_pos,n_clust) = squeeze(iage(lon_pos,edge(i,lon_pos),i));
        peak_period_vec(i,data_pos,n_clust) = squeeze(peak_period(lon_pos,edge(i,lon_pos),i));
        fsdrad_vec(i,data_pos,n_clust) = squeeze(fsdrad(lon_pos,edge(i,lon_pos),i));

        % ADD evolution of the FSD Terms
        dafsd_latg_vec1(i,data_pos,n_clust) = squeeze(dafsd_latg(lon_pos,edge(i,lon_pos),1,i));
        dafsd_latm_vec1(i,data_pos,n_clust) = squeeze(dafsd_latm(lon_pos,edge(i,lon_pos),1,i));
        dafsd_newi_vec1(i,data_pos,n_clust) = squeeze(dafsd_newi(lon_pos,edge(i,lon_pos),1,i));
        dafsd_weld_vec1(i,data_pos,n_clust) = squeeze(dafsd_weld(lon_pos,edge(i,lon_pos),1,i));
        dafsd_wave_vec1(i,data_pos,n_clust) = squeeze(dafsd_wave(lon_pos,edge(i,lon_pos),1,i));

    end
    
end
end







%% plotting ----------------------------------------------------------------
% Data must be timeseries
cluster_idx = 2;

data.lat = mean(lat_vec(:,:,cluster_idx),2);
data.aice = mean(aice_vec(:,:,cluster_idx),2);
data.windspeed = mean(windspeed_vec(:,:,cluster_idx),2);
data.icespeed = mean(icespeed_vec(:,:,cluster_idx),2);
data.winddir = mean(winddir_vec(:,:,cluster_idx),2);
data.icedir = mean(icedir_vec(:,:,cluster_idx),2);
data.Tair = mean(Tair_vec(:,:,cluster_idx),2);
data.sst = mean(sst_vec(:,:,cluster_idx),2);
data.peak_period = mean(peak_period_vec(:,:,cluster_idx),2);
data.swh = mean(swh_vec(:,:,cluster_idx),2);
data.daidtd = mean(daidtd_vec(:,:,cluster_idx),2);
data.daidtt = mean(daidtt_vec(:,:,cluster_idx),2);
data.fsdrad = mean(fsdrad_vec(:,:,cluster_idx),2);
data_dafsd(1,:) = mean(dafsd_latg_vec1(:,:,cluster_idx),2);
data_dafsd(2,:) = mean(dafsd_latm_vec1(:,:,cluster_idx),2);
data_dafsd(3,:) = mean(dafsd_newi_vec1(:,:,cluster_idx),2);
data_dafsd(4,:) = mean(dafsd_weld_vec1(:,:,cluster_idx),2);
data_dafsd(5,:) = mean(dafsd_wave_vec1(:,:,cluster_idx),2);
%               mean(dafsd_latm_vec1(:,:,cluster_idx),2),
%               mean(dafsd_newi_vec1(:,:,cluster_idx),2),
%               mean(dafsd_weld_vec1(:,:,cluster_idx),2),
%               mean(dafsd_wave_vec1(:,:,cluster_idx),2)];

n_clust_vec = 2;
plotting = "on";
if plotting == "on"
    close all
    clc
    conFigure(11,4)
    figure
    %sgtitle(sprintf('Lon %g E',lon(lon_pos,1)))
    sgtitle(sprintf('Mean values for Class %g',clust_vec(cluster_idx)))
end

if plotting == "on"
    subplot(8,1,1)
    title(sprintf('Class %g',clust_vec(cluster_idx)))
    yyaxis left
    stairs(date_label(1:end-1),lat(lon_pos,edge(:,10)))
    ylabel('Edge location [$^\circ$S]')
    ylim([-70,-50])
    yyaxis right
    plot(data.aice)
    ylim([0,1.2])
    ylabel('SIC [$\%$]')
    
    subplot(8,1,2)
    yyaxis left
    plot(date_label(1:end-1),data.windspeed)
    ylabel('Wind speed [m/s]')
    %ylim([-30,20])
    yyaxis right
    plot(data.icespeed)
    ylabel('Ice speed [m/s]') % atan(0) = N
    %ylim([0,1])
    %hold off
%     
%     north_wind = abs(winddir_vec) < pi/4;
%     south_wind = abs(winddir_vec) > pi/4;
%     bi_wind = north_wind - south_wind;
% 
%     subplot(5,num_clusters,9+n_clust-1)
%     yyaxis left
%     stairs(delta_edge(:,lon_pos),'LineWidth',1,'LineStyle','--')
%     ylim([-3,3])
%     ylabel('$\Delta$ edge')
%     yyaxis right
%     stairs(bi_wind,'LineWidth',1,'LineStyle','--')
%     ylabel('N/S Wind')
%     ylim([-3,3])
    %corr(delta_edge(:,lon_pos), bi_wind(1:end-1)')
    %edge_idx_vec = edge_idx(:,lon_pos);
    %same = delta_edge(:,lon_pos).*edge_idx_vec == bi_wind(1:end-1)'.*edge_idx_vec;
    %percent_wind = sum(same)/length(delta_edge(:,lon_pos));
    %sprintf('%g percent of the time wind direction correlates with change in edge location',round(percent_wind.*100))

    subplot(8,1,3)
    yyaxis left
    plot(date_label(1:end-1),data.winddir.*180)
    ylabel('Wind direction [$^\circ$]')
    %ylim([-200,200+360])
    yyaxis right
    plot(data.icedir.*180)
    ylabel('Ice direction [$^\circ$]') % atan(0) = N
    %ylim([-200-360,200])


    subplot(8,1,4)
    yyaxis left
    plot(date_label(1:end-1),data.Tair)
    ylabel('Air temp. [$^\circ$C]')
    ylim([-15,5])
    yyaxis right
    plot(data.sst)
    ylabel('SST [$^\circ$C]') % atan(0) = N
    ylim([-2,0])

    north_ice = abs(icedir_vec) < pi/4;
    south_ice = abs(icedir_vec) > pi/4;
    bi_ice = north_ice - south_ice;

    subplot(8,1,5)
    yyaxis left
    %stairs(delta_edge(:,lon_pos),'LineWidth',1,'LineStyle','--')
    %ylim([-3,3])
    %ylabel('$\Delta$ edge')
    plot(date_label(1:end-1),data.peak_period)
    ylim([0,20])
    ylabel('Peak period [s]')
    yyaxis right
    plot(data.swh)
    ylim([0,9])
    ylabel('SWH [m]')

    subplot(8,1,6)
    yyaxis left
    plot(date_label(1:end-1),data.daidtd)
    yline(0)
    ylabel('Daidtd [$\%$]') 
    ylim([-100,100])
    yyaxis right
    plot(data.daidtt)
    ylabel('Daidtt [$\%$]') 
    ylim([-100,100])

    subplot(8,1,7)
    yyaxis left
    plot(date_label(1:end-1),data.fsdrad)
    ylabel('Mean FSD [m]') 
    ylim([0,100])
    %yyaxis right

    subplot(8,1,8)
    yyaxis left
    plot(data_dafsd(3,:)')
    ylabel('Change in FSD1 [$\%$]') 
    ylim([-10^(-2),10^(-2)])
    yline(0)
    yyaxis right
    plot(data_dafsd([1,2,4,5],:)')
    ylim([-10^(-3),10^(-3)])
    legend("Newi","","Latg","Latm","Weld","Wave")

    %ylim([-100,100])
    
    %corr_direction = [];
    %corr_direction(n_clust) = corr(winddir_vec',icedir_vec');
end

%% Calculate the changes in ice edge in the wind direction
close all
clc
[lat,lon] = grid_read('om2');
for lon_pos = 1:length(edge_idx(1,:))
for n_clust = 1:num_clusters
    for i = 1:n_files-1
        [lat_ice_edge_new, lon_ice_edge_new, edge(i,:)] = find_edge(k_means_cube(:,:,i),n_clust,sector,lat,lon);
    end
    [~,wid] = size(edge);
    for i = 1:n_files-2
        % Difference in ice edge
        delta_edge(i,:) = edge(i+1,:) - edge(i,:);
        for j = 1:length(delta_edge(i,:))
            if delta_edge(i,j) > 0
                edge_idx(i,j) = 1;
            end
        end
    end
    swh_vec = [];
    winddir_vec = [];
    for i = 1:length(edge(:,lon_pos)) % Each file
        swh_vec(i) = squeeze(swh(lon_pos,edge(i,lon_pos),i));
        windspeed_vec(i) = sqrt(squeeze(uatm(lon_pos,edge(i,lon_pos),i)).^2 + squeeze(vatm(lon_pos,edge(i,lon_pos),i)).^2);
        winddir_vec(i) = atan(squeeze(vatm(lon_pos,edge(i,lon_pos),i))/squeeze(uatm(lon_pos,edge(i,lon_pos),i)));
    end
    
    edge_idx_vec = edge_idx(:,lon_pos);
    same = delta_edge(:,lon_pos).*edge_idx_vec == bi_wind(1:end-1)'.*edge_idx_vec;
    percent_wind(n_clust,lon_pos) = sum(same)/length(delta_edge(:,lon_pos));
    %sprintf('%g percent of the time wind direction correlates with change in edge location',round(percent_wind.*100))
end
end

%% Correlation of wind direction at these times
region_label = {'Large','Consolidated','Low SIC','Pancake'};
close all
conFigure(11,3)
figure
plot(lon(1:length(edge_idx(1,:)),1),percent_wind)
xlabel('Longitude [$^\circ$E]')
ylabel('Agreement between wind direction and change of edge location')
ylim([0,1])
legend(region_label)



%% FUNCTIONS

function [lat_ice_edge_new, lon_ice_edge_new, edge] = find_edge(data,class,sector,lat,lon)
%FIND_ICE_EDGE finds the ice edge of ice area data aice according to a
%threshold SIC
%   Detailed explanation goes here
    [len,wid] = size(data);
    ice_edge = zeros(len,wid);
    
    coords = sector_coords(sector);
    [lat_north, lon_east] = lat_lon_finder(coords(1,1),coords(1,2),lat,lon);
    [lat_south, lon_west] = lat_lon_finder(-75,coords(3,2),lat,lon);

    edge = ones(1,len);
    
    if lon_east > lon_west
        for j = lon_east:len
            for i = lat_north:-1:lat_south
                if data(j,i) == class
                    edge(j) = i+1;
                    break
                end
            end
        end
        for j = 1:lon_west
            for i = lat_north:-1:lat_south
                if data(j,i) == class
                    edge(j) = i+1;
                    break
                end
            end
        end

    else
        for j = lon_east:lon_west
            for i = lat_north:-1:lat_south
                if data(j,i) == class
                    edge(j) = i+1;
                    break
                end
            end
        end
    end
    
    for j = 1:len
        lat_ice_edge(j) = lat(j,edge(j));
        lon_ice_edge(j) = lon(j,edge(j));
    end
    
    % Make the ice edge fit the resolution
    count = 2;
    lat_ice_edge_new(1) = lat_ice_edge(1);
    lat_ice_edge_new(1) = lat_ice_edge(1);
    for j = 2:len
        if lat_ice_edge(j) ~= lat_ice_edge(j-1)
            % Then we have a step
            lat_ice_edge_new(count) = lat_ice_edge(j);
            lon_ice_edge_new(count) = lon_ice_edge(j-1);
            count = count + 1;
        end
        lat_ice_edge_new(count) = lat_ice_edge(j);
        lon_ice_edge_new(count) = lon_ice_edge(j);
        count = count + 1;
    end

end