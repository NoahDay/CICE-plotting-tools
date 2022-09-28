% Cyclones were manually identified using a threshold of a minimum isobar
% value of 940 for 2019. 4 cyclones/storms were identified that meet these
% requirements. From which, we will calculate the change in the ice edge
% and consolidated boundaries to determine whether these correlate with
% firstly wind direction and ice velocities. Next, we will investigate the
% role of thermodynamics, i.e., 
% Did the surface temperature increase? 
% Is the cyclone pushing hot air further south and bringing colder air north? 
% Is new ice forming and where?
% Is wave height increasing over these times, and are the waves breaking up
% ice? and or promoting the formation of pancakes rather than nilas ice?


% Add libraries
clear all
addpath functions
sector = "EA";
% Setup
historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/2019/may-sep/cyclone940hpa_4/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end

for i = 1:365
    temp = datetime(dirdates(1,:)) - i;
    if temp == "01-May-2019"
        date_diff = i;
    end
end
%date_diff = 106; % There are 106 days between the start of May and the cyclone hitting
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
%
[lat,lon] = grid_read('om2');
for file_number = date_diff:date_diff+n_files-1
    row_file = find(index_both);
    row_vec = row_file(file_number)+1:row_file(file_number+1);
    file_idx = idx(row_vec);
    [~, sector_mask] = data_format_sector(filenames(file_number-date_diff+1,:),'aice',sector);
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

for i = 1:n_files
    k_means_cyclone(:,:,i) = k_means_cube(:,:,date_diff+i-1);
end
clear k_means_cube index_lat index_lon index_both file_idx row_vec row_idx k_means X_new X_map kmeans_cluster label_vec
%% Get boundary data
clear swh_vec windspeed_vec winddir_vec icespeed_vec icedir_vec Tair_vec sst_vec daidtt_vec daidtd_vec aice_vec iage_vec peak_period_vec edge lat_vec
date_label = datetime(dirdates);
%lon_pos = 10;
clust_vec = [2,4,5];
n_clust_vec = length(clust_vec);

longs = 10:140;
for lon_pos = longs
for n_clust = 1:n_clust_vec
    for i = 1:n_files
        [lat_ice_edge_new, lon_ice_edge_new, edge(i,:)] = find_edge(k_means_cyclone(:,:,i),clust_vec(n_clust),sector,lat,lon);
        [lat_ice_edge_new, lon_ice_edge_new, ice_edge(i,:)] = find_ice_edge(aice(:,:,i),0.15,sector,lat,lon);
        if clust_vec(n_clust) == 5
            edge(i,:) = ice_edge(i,:)-1;
            
        end
        edge_stored(i,:,n_clust) = edge(i,:);
        %if n_clust == 2
        %edge(i,:) = edge(i,:) - 1;
        %end
    end

    [~,wid] = size(edge);
    for i = 1:n_files-1
        % Difference in ice edge
        delta_edge(i,:,n_clust) = edge(i+1,:) - edge(i,:);
        delta_ice_edge(i,:) = ice_edge(i+1,:) - ice_edge(i,:);
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
        vvel_vec(i,data_pos,n_clust) =squeeze(vvel(lon_pos,edge(i,lon_pos),i));
        vatm_vec(i,data_pos,n_clust) =squeeze(vatm(lon_pos,edge(i,lon_pos),i));
        Tair_vec(i,data_pos,n_clust) = squeeze(Tair(lon_pos,edge(i,lon_pos),i));
        sst_vec(i,data_pos,n_clust) = squeeze(sst(lon_pos,edge(i,lon_pos),i));
        daidtt_vec(i,data_pos,n_clust) = squeeze(daidtt(lon_pos,edge(i,lon_pos),i));
        daidtd_vec(i,data_pos,n_clust) = squeeze(daidtd(lon_pos,edge(i,lon_pos),i));
        frzmlt_vec(i,data_pos,n_clust) = squeeze(frzmlt(lon_pos,edge(i,lon_pos),i));
        
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

%% Plot the change in the boundaries
spacing = 0.2;
[X,Y] = meshgrid(-2:spacing:2);
Z = X.*exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,spacing);
clf
figure
quiver(X,Y,DX,DY)
% Add quivers for wind direction

%%
lineWidth = 1.5;
xrange = [10,150];
close all
%clf

conFigure(11)
figure('Position',[1,1,20,20])
sgtitle('Change in consolidated boundary')
for i = 1:n_files-1
    subplot(n_files-1,4,4*i-3)
    yyaxis left
    stairs(lon(1:end,1),delta_edge(i,:,1),'LineWidth',lineWidth)
    title(dirdates(i+1,:))
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[vatm_vec(i+1,:,1); vvel_vec(i+1,:,1)*50],'LineWidth',lineWidth/2)
    ylim([-30,30])
    ylabel('[ms$^{-1}$]')
    if i == 1
        legend("","",'$v_{wind}$','$50\times v_{ice}$','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')

    idx_frz = frzmlt_vec(i+1,:,1) > 0;
    frz_vec = zeros(1,length(frzmlt_vec(i+1,:,1)));
    frz_vec(idx_frz) = frzmlt_vec(i+1,idx_frz,1);

    subplot(n_files-1,4,4*i-2)
    yyaxis left
    stairs(lon(1:end,1),delta_edge(i,:,1),'LineWidth',lineWidth)
    title(dirdates(i+1,:))
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[Tair_vec(i+1,:,1); sst_vec(i+1,:,1); frz_vec],'LineWidth',lineWidth/2)
    ylim([-30,30])
    ylabel('[$^\circ$C]')
    
    if i == 1
        legend("","",'$T_{air}$','SST','Frz pot.','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')


    subplot(n_files-1,4,4*i-1)
    yyaxis left
    stairs(lon(1:end,1),delta_edge(i,:,1),'LineWidth',lineWidth)
    title(dirdates(i+1,:))
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[daidtt_vec(i+1,:,1); daidtd_vec(i+1,:,1)],'LineWidth',lineWidth/2)
    ylim([-100,100])
    ylabel('[SIC $\%$]')
    if i == 1
        legend("","",'Thermo','Dynamics','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')

    subplot(n_files-1,4,4*i)
    title(dirdates(i+1,:))
    ylim([-0.4,0.4])
    xlim(xrange)
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    hold on
    plot(lon(longs,1),[dafsd_latg_vec1(i+1,:,1); dafsd_latm_vec1(i+1,:,1); dafsd_newi_vec1(i+1,:,1); dafsd_weld_vec1(i+1,:,1); dafsd_wave_vec1(i+1,:,1);],'LineWidth',lineWidth)
    hold off
    grid on
    %ylim([-100,100])
    ylabel('dafsd1 [$\%$]')
    if i == 1
        lgd = legend("","Latg","Latm","Newi","Weld","Wave",'Location','south','Orientation','horizontal');
        lgd.NumColumns = 3;
    end
    xlabel('Latitude [$^\circ$E]')

     
end


%%
%close all
figure('Position',[10,1,20,20])
conFigure(11)
for i = 1:n_files-1
    subplot(n_files-1,4,4*i-3)
    sgtitle('Change in ice edge boundary')
    yyaxis left
    stairs(lon(1:end,1),delta_ice_edge(i,:),'LineWidth',lineWidth)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    title(dirdates(i+1,:))
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[vatm_vec(i+1,:,3); vvel_vec(i+1,:,3)*50])
    ylim([-30,30])
    ylabel('[ms$^{-1}$]')
    if i == 1
        legend("","",'$v_{wind}$','$50\times v_{ice}$','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')

    idx_frz = frzmlt_vec(i+1,:,3) > 0;
    frz_vec = zeros(1,length(frzmlt_vec(i+1,:,3)));
    frz_vec(idx_frz) = frzmlt_vec(i+1,idx_frz,3);

    subplot(n_files-1,4,4*i-2)
    yyaxis left
    stairs(lon(1:end,1),delta_ice_edge(i,:,1),'LineWidth',lineWidth)
    title(dirdates(i+1,:))
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[Tair_vec(i+1,:,3); sst_vec(i+1,:,3); frz_vec],'LineWidth',lineWidth/2)
    ylim([-30,30])
    ylabel('[$^\circ$C]')
    if i == 1
        legend("","",'$T_{air}$','SST','Frz Pot.','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')

    subplot(n_files-1,4,4*i-1)
    yyaxis left
    stairs(lon(1:end,1),delta_ice_edge(i,:,1),'LineWidth',lineWidth)
    title(dirdates(i+1,:))
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    ylabel('Change in boundary')
    grid on
    yticks(-5:5)
    ylim([-5,5])
    xlim(xrange)
    yyaxis right
    plot(lon(longs,1),[daidtt_vec(i+1,:,3); daidtd_vec(i+1,:,3)],'LineWidth',lineWidth/2)
    ylim([-100,100])
    ylabel('[SIC $\%$]')
    if i == 1
        legend("","",'Thermo','Dynamics','Location','south','Orientation','horizontal')
    end
    yyaxis left
    xlabel('Latitude [$^\circ$E]')

    subplot(n_files-1,4,4*i)
    title(dirdates(i+1,:)) 
    ylim([-0.4,0.4])
    xlim(xrange)
    yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
    hold on
    plot(lon(longs,1),[dafsd_latg_vec1(i+1,:,3); dafsd_latm_vec1(i+1,:,3); dafsd_newi_vec1(i+1,:,3); dafsd_weld_vec1(i+1,:,3); dafsd_wave_vec1(i+1,:,3);],'LineWidth',lineWidth)
    hold off
    grid on
    %ylim([-100,100])
    ylabel('dafsd1 [$\%$]')
    if i == 1
        lgd = legend("","Latg","Latm","Newi","Weld","Wave",'Location','south','Orientation','horizontal');
        lgd.NumColumns = 3;
    end
    xlabel('Latitude [$^\circ$E]')
end

% subplot(3,1,2)
% stairs(lon(1:end,1),delta_ice_edge(2,:),'LineWidth',lineWidth)
% yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
% title('Day 2')
% ylim([-5,5])
% xlim(xrange)
% yyaxis right
% plot(lon(longs,1),[vatm_vec(3,:,3); vvel_vec(3,:,3)*50])
% ylim([-30,30])
% yyaxis left
% 
% 
% subplot(3,1,3)
% stairs(lon(1:end,1),delta_ice_edge(3,:),'LineWidth',lineWidth)
% yline(0,'LineWidth',lineWidth,'LineStyle','--','color','k','Alpha',0.5)
% title('Day 3')
% ylim([-5,5])
% xlim(xrange)
% xlabel('Longitude')
% yyaxis right
% plot(lon(longs,1),[vatm_vec(4,:,3); vvel_vec(4,:,3)*50])
% %plot(lon(longs,1),aice_vec(4,:,3))
% ylim([-30,30])
% yyaxis left


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
                    edge(j) = i;
                    break
                end
            end
        end
        for j = 1:lon_west
            for i = lat_north:-1:lat_south
                if data(j,i) == class
                    edge(j) = i;
                    break
                end
            end
        end

    else
        for j = lon_east:lon_west
            for i = lat_north:-1:lat_south
                if data(j,i) == class
                    edge(j) = i;
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
