%% Single cell analysis
addpath functions
clear all 
clc
historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/2019/';%'/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history_h/';
sector = "SH";
user = "noahday";
grid = 'om2';
%% Read in the data
a = dir([historydir '/*.nc']);
n_files = numel(a);
[lat,lon,~,ulat,ulon] = grid_read('om2');
[Lat,Lon] = meshgrid(lat(1,:),lon(:,180));
SIC = 0.15;
% Read in file names
for i = 1:n_files
    filename.cice(i,:) = strcat(historydir,a(i).name);
    dirdates(i,:) = a(i).name(11:end-3);
    if i == 1
        NFSD = ncread(filename.cice(1,:),"NFSD");
        [floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
    end

    [aice(:,:,i), sector_mask] = data_format_sector(filename.cice(i,:), "aice",sector);
    [lat_ice_edge, lon_ice_edge, edge(i,:)] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
    air_u(:,:,i) = data_format_sector(filename.cice(i,:),"uatm",sector); % m/s
    air_v(:,:,i) = data_format_sector(filename.cice(i,:),"vatm",sector); % m/s
    ice_u(:,:,i) = data_format_sector(filename.cice(i,:),"uvel",sector); % m/s
    ice_v(:,:,i) = data_format_sector(filename.cice(i,:),"vvel",sector); % m/s
    daidtd(:,:,i) = data_format_sector(filename.cice(i,:),"daidtd",sector); % SIC
    daidtt(:,:,i) = data_format_sector(filename.cice(i,:),"daidtt",sector); % SIC
    frzmlt(:,:,i) = data_format_sector(filename.cice(i,:),"frzmlt",sector); % 

    wave_sig_ht(:,:,i) = data_format_sector(filename.cice(i,:),"wave_sig_ht",sector); % m
    fsdrad(:,:,i) = data_format_sector(filename.cice(i,:),"fsdrad",sector); % m
    sst(:,:,i) = data_format_sector(filename.cice(i,:),"sst",sector); % C
    airtmp(:,:,i) = data_format_sector(filename.cice(i,:),"Tair",sector); % C


    afsdn(:,:,:,:) = data_format_sector(filename.cice(i,:),"afsdn",sector); % m
    newi(:,:,:) = data_format_sector(filename.cice(i,:),"dafsd_newi",sector); % m
    weld(:,:,:) = data_format_sector(filename.cice(i,:),"dafsd_weld",sector); % m
    latg(:,:,:) = data_format_sector(filename.cice(i,:),"dafsd_latg",sector); % m
    latm(:,:,:) = data_format_sector(filename.cice(i,:),"dafsd_latm",sector); % m
    wave(:,:,:) = data_format_sector(filename.cice(i,:),"dafsd_wave",sector); % m
    for j = 1:numel(NFSD)
        itd1(:,:,j,i) = afsdn(:,:,j,1).*floe_binwidth(j);

        dafsd_newi(:,:,j,i) = newi(:,:,j).*floe_binwidth(j);
        dafsd_weld(:,:,j,i) = weld(:,:,j).*floe_binwidth(j);
        dafsd_latg(:,:,j,i) = latg(:,:,j).*floe_binwidth(j);
        dafsd_latm(:,:,j,i) = latm(:,:,j).*floe_binwidth(j);
        dafsd_wave(:,:,j,i) = wave(:,:,j).*floe_binwidth(j);
    end

    [Lat,Lon,u,v] = recenter(Lat,Lon,ice_u(:,:,i),ice_v(:,:,i));
    Curl(:,:,i) = cdtcurl(Lat,Lon,u,v);
    % Multiply in the south with -1
    Curl(:,:,i) = Curl(:,:,i).*sign(Lat);

    [Lat,Lon,u,v] = recenter(Lat,Lon,air_u(:,:,i),air_v(:,:,i));
    air_curl(:,:,i) = cdtcurl(Lat,Lon,u,v);
    % Multiply in the south with -1
    air_curl(:,:,i) = air_curl(:,:,i).*sign(Lat);
    if mod(i,10) == 1
        disp(i)
    end
end


% Grid type from COSIMA
hte = ncread(filename.cice(i,:),'HTE');
row = 281;
% OM2 grid
hte = rearrange_matrix(hte,row,2);
hte = [hte(1,:); hte];
%% Air pressure
filename.airPressure = '/Users/noahday/GitHub/CICE-plotting-tools/observations/jra55/psl_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-4-0_gr_201701010000-201712312100.nc';
airPressure = ncread(filename.airPressure, "psl");
lat_jra = ncread(filename.airPressure, "lat");
lon_jra = ncread(filename.airPressure, "lon");
sector = "AU";
coords = sector_coords(sector);
% September start iceh_inst.2007-07-01-75600.nc
jra_idx = 120*8; % May 1st
%181*8; % July 1st
[len,wid,~] = size(airPressure);
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
                    %airPressureInterp(j,k,1:n_files) = interp1(1:ceil(n_files/3),temp,linspace(1,ceil(n_files/3),n_files));
                    airPressureInterp(j,k,i) = temp;
                else
                    airPressureInterp(j,k,i) = NaN;
                end
            end
        end
    end
end
%% Plot ice and air velocities for a cell
% (lat = -62, lon = 8)
% FIX Ice edge
lat_desired = -61;
lon_desired = 2;

lat_pos_jra = near1(lat_jra,lat_desired);
lon_pos_jra = near1(lon_jra, lon_desired);


[lon_pos,lat_pos,~] = near2(lon,lat,lon_desired,lat_desired);

%lat_pos = edge(:,lon_pos);
centre_cyclone = find(min(squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra,:))) == squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra,:)));

pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
temp = aice(:,:,1);
temp(lon_pos,lat_pos(2)) = 10;
font_size = 12;
close all
figure
w = worldmap('world');
    axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [0 0 0]);
    setm(w, 'maplatlimit', [-80,-30]);
    setm(w, 'maplonlimit', [345,90]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,temp)
    ylabel('Latitude')
    scalebar('color',[0.0 0.0 0.0], 'location','sw');
    cb = colorbar; set(cb,'position',[.89 .15 .02 .75])
    cb.Label.String = 'CURL '; cb.Label.Interpreter = 'latex';
    cb.FontSize = font_size;cb.Color = 'black';
    cmocean('balance',31)
    caxis([0,15])
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    %[c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LevelStep',20);%,'ShowText','on')
    %quivermc(ulat,ulon,ice_u(:,:,i),ice_v(:,:,i),'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
    %t = clabelm(c1,h);
    %set(t,'Color','r'); set(t,'BackgroundColor','none'); 
    %set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])

% Ice speed
ice_mag_vec(:) = sqrt(ice_u(lon_pos,lat_pos,:).^2 + ice_v(lon_pos,lat_pos,:).^2);
air_mag_vec(:) = sqrt(air_u(lon_pos,lat_pos,:).^2 + air_v(lon_pos,lat_pos,:).^2);



%% Plotting
close all
n_panes = 4;
conFigure(11,0.5)
figure
subplot(n_panes,2,1)
yyaxis left
plot(ice_mag_vec)
ylabel('Ice speed [m/s]')
hold on
yyaxis right
plot(air_mag_vec)
ylabel('Wind speed [m/s]')
xline(centre_cyclone,'LineStyle','--')
hold off


% Lat_pos = near1(Lat(1,:),-62);
% Lon_pos = near1(Lon(:,180), 8);
% curl_vec(:) = Curl(Lon_pos,Lat_pos,:);
% air_curl_vec(:) = air_curl(Lon_pos,Lat_pos,:);
% figure
% yyaxis left
% plot(curl_vec)
% ylabel('curl of ice')
% hold on
% yyaxis right
% plot(air_curl_vec)
% ylabel('Curl of wind')
% hold off


% Change in SIC
daidtd_vec(:) = daidtd(lon_pos,lat_pos,:);
daidtt_vec(:) = daidtt(lon_pos,lat_pos,:);
SIC_vec = squeeze(aice(lon_pos,lat_pos,:));
dt = 3600;
dSIC_vec = [0; SIC_vec(2:end)-SIC_vec(1:end-1)].*100;
% dSIC_cice = daidtd_vec + daidtt_vec;
% 
% daidtt_noah(1) = SIC_vec(1);
% for i = 2:numel(SIC_vec)
%     daidtt_noah(i) = SIC_vec(i) - daidtt_noah(i-1)/dt;
% end
% daidtt_noah = daidtt_noah.*100;
% 
% RMSE = sqrt(mean((dSIC_cice - daidtt_noah).^2))

%idx = daidtt_vec > 100;
%daidtt_vec(idx == 100);
subplot(n_panes,2,3)
%yyaxis left
plot(daidtd_vec)
ylabel('Change in SIC [SIC/model step]')
ylim([-100,100])

hold on
plot(daidtt_vec)
plot(dSIC_vec)
xline(centre_cyclone,'LineStyle','--')
legend(["Dynamics","Thermo.","Total change",""],'Location','northwest')
%yyaxis right

%ylabel('Change in SIC from thermodynamics')
%ylim([-100,100])

hold off



% Distance from ice edge and aice
%south_point = find(MIZ_vec_temp, 1, 'first');
%north_point = find(MIZ_vec_temp, 1, 'last');
%MIZ_width(i,:) = pathdist([ulat_vec_temp(south_point-1),ulat_vec_temp(north_point)],[ulon_vec_temp(south_point),ulon_vec_temp(north_point)],'km');
for i = 1:n_files
    dist_from_edge(i) = sum(hte(lon_pos,lat_pos:edge(i,lon_pos)));
end

subplot(n_panes,2,5)
yyaxis left
plot(squeeze(aice(lon_pos,lat_pos,:)))
ylabel('SIC')
ylim([0,1.2])
hold on
yyaxis right
plot(dist_from_edge/1000)
ylabel('Distance from ice edge [km]')
%ylim([0,500])
xline(centre_cyclone,'LineStyle','--')
hold off


subplot(n_panes,2,2)
yyaxis left
plot(squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra,:))/100)
ylabel('Air pressure [hPa]')
%ylim([0,1.2])
hold on
yyaxis right
plot(squeeze(airtmp(lon_pos,lat_pos,:)))
ylabel('Air temperature [C]')
%ylim([0,4])
xline(centre_cyclone,'LineStyle','--')
hold off


subplot(n_panes,2,4)
yyaxis left
plot(squeeze(sst(lon_pos,lat_pos,:)))
ylabel('Sea surface temperature [C]')
%ylim([0,1.2])
hold on
yyaxis right
plot(squeeze(frzmlt(lon_pos,lat_pos,:)))
ylabel('Freeze/melt potential')
%ylim([0,4])
xline(centre_cyclone,'LineStyle','--')
hold off


data = [squeeze(dafsd_latg(lon_pos,lat_pos,1,:)), squeeze(dafsd_latm(lon_pos,lat_pos,1,:)), squeeze(dafsd_newi(lon_pos,lat_pos,1,:)), squeeze(dafsd_weld(lon_pos,lat_pos,1,:)), squeeze(dafsd_wave(lon_pos,lat_pos,1,:))];

subplot(n_panes,2,6)
%yyaxis left
area(data)
ylabel('Change in FSD cat 1 [\%]')

%ylim([0,1.2])
% hold on
% plot(squeeze(dafsd_latm(lon_pos,lat_pos,1,:)))
% plot(squeeze(dafsd_newi(lon_pos,lat_pos,1,:)))
% plot(squeeze(dafsd_weld(lon_pos,lat_pos,1,:)))
% plot(squeeze(dafsd_wave(lon_pos,lat_pos,1,:)))
% % yyaxis right
% plot(squeeze(pancake_weld(lon_pos,lat_pos,:)))
% ylabel('Pancake welding [\%]')
%ylim([0,1.2])
xline(centre_cyclone,'LineStyle','--')
legend(["Latg","Latm","Newi","Weld","Wave",""],'Location','northwest')
% hold off


subplot(n_panes,1,4)
%yyaxis left
area(squeeze(itd1(lon_pos,lat_pos,:,:))')
ylabel('FSD SIC [\%]')

%ylim([0,1.2])
%yyaxis right
%plot(squeeze(pancake(lon_pos,lat_pos,:)))
%ylabel('')
%ylim([0,1.2])
xline(centre_cyclone,'LineStyle','--')
legend([string(1:12),""],'Location','northwest','Orientation','horizontal','NumColumns',3)


% subplot(n_panes,1,7)
% yyaxis left
% plot(squeeze(wave_sig_ht(lon_pos,lat_pos,:)))
% ylabel('SWH [m]')
% %ylim([0,1.2])
% hold on
% yyaxis right
% plot(squeeze(pancake_weld(lon_pos,lat_pos,:)))
% ylabel('Pancake welding [\%]')
% %ylim([0,1.2])
% xline(centre_cyclone,'LineStyle','--')
% hold off


% SST
% Air temp
% Tsfc
%% Plotting ice edge
lat_desired = -61;
lon_desired = 2;

lat_pos_jra = near1(lat_jra,lat_desired);
lon_pos_jra = near1(lon_jra, lon_desired);
lat_diff = -1;

[lon_pos,lat_pos,~] = near2(lon,lat,lon_desired,lat_desired);
lat_pos = edge(:,lon_pos) + lat_diff;
centre_cyclone = find(min(squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra,:))) == squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra,:)));
centre_cyclone = 0;
dt = 3600; % CICE timestep

% Read in data along ice edge
for i = 1:n_files
    lat_pos_jra(i) = near1(lat_jra, lat(lon_pos,edge(i,lon_pos))) + lat_diff;
    daidtd_vec(i) = squeeze(daidtd(lon_pos,lat_pos(i),i));
    daidtt_vec(i) = daidtt(lon_pos,lat_pos(i),i);
    SIC_vec(i) = squeeze(aice(lon_pos,lat_pos(i),i));
    pressure_vec(i) = squeeze(airPressureInterp(lon_pos_jra,lat_pos_jra(i),i))/100;

    ice_mag_vec(i) = sqrt(ice_u(lon_pos,lat_pos(i),i).^2 + ice_v(lon_pos,lat_pos(i),i).^2);
    air_mag_vec(i) = sqrt(air_u(lon_pos,lat_pos(i),i).^2 + air_v(lon_pos,lat_pos(i),i).^2);

    dist_from_edge(i) = sum(hte(lon_pos,lat_pos(i):edge(i,lon_pos)));
    airtmp_vec(i) = squeeze(airtmp(lon_pos,lat_pos(i),i));
    sst_vec(i) = squeeze(sst(lon_pos,lat_pos(i),i));
    frzmlt_vec(i) = squeeze(frzmlt(lon_pos,lat_pos(i),i));

    dafsd_data(i,:) = [squeeze(dafsd_latg(lon_pos,lat_pos(i),1,i)), squeeze(dafsd_latm(lon_pos,lat_pos(i),1,i)), squeeze(dafsd_newi(lon_pos,lat_pos(i),1,i)), squeeze(dafsd_weld(lon_pos,lat_pos(i),1,i)), squeeze(dafsd_wave(lon_pos,lat_pos(i),1,i))];
    itd1_data(i,:) = squeeze(itd1(lon_pos,lat_pos(i),:,i))';
    lat_vec(i) = lat(lon_pos,lat_pos(i));
    if i == 1
        xlabel_date(i) = datetime(2019,5,1);
    else
        xlabel_date(i) = xlabel_date(i-1)+1;
    end
end
dSIC_vec = [0, SIC_vec(2:end)-SIC_vec(1:end-1)].*100;
delta_lat = [0, lat_vec(2:end) - lat_vec(1:end-1)];
lat_idx = delta_lat ~=0;
SIC_change = lat_idx.*(1:numel(lat_idx));
SIC_change(1) = 0;
SIC_change = SIC_change(lat_idx)-1;
close all
n_panes = 4;
conFigure(11,0.5)
figure
subplot(n_panes,2,1)
yyaxis left
plot(xlabel_date,ice_mag_vec)
ylabel('Ice speed [m/s]')
hold on
yyaxis right
plot(xlabel_date,air_mag_vec)
ylabel('Wind speed [m/s]')
xline(SIC_change,'LineStyle','--')

hold off



subplot(n_panes,2,3)
%yyaxis left
plot(xlabel_date,daidtd_vec)
ylabel('Change in SIC [SIC/model step]')
ylim([-100,100])
hold on
plot(xlabel_date,daidtt_vec)
%plot(xlabel_date,dSIC_vec)
xline(SIC_change,'LineStyle','--')

legend(["Dynamics","Thermo.","",""],'Location','northwest')
%yyaxis right

%ylabel('Change in SIC from thermodynamics')
%ylim([-100,100])

hold off



% Distance from ice edge and aice

subplot(n_panes,2,5)
yyaxis left
plot(xlabel_date,SIC_vec)
%yline(0.15)
ylabel('SIC')
ylim([0,1])
hold on
yyaxis right
plot(xlabel_date,lat_vec)
ylabel('Latitude [$^\circ$]')
%ylim([0,500])
xline(SIC_change,'LineStyle','--')
hold off


subplot(n_panes,2,2)
yyaxis left
plot(xlabel_date,pressure_vec)
ylabel('Air pressure [hPa]')
%ylim([0,1.2])
hold on
yyaxis right
plot(xlabel_date,airtmp_vec)
ylabel('Air temperature [C]')
%ylim([0,4])
xline(SIC_change,'LineStyle','--')
hold off


subplot(n_panes,2,4)
yyaxis left
plot(xlabel_date,sst_vec)
ylabel('Sea surface temperature [C]')
%ylim([0,1.2])
hold on
yyaxis right
plot(xlabel_date,frzmlt_vec)
ylabel('Freeze/melt potential')
%ylim([0,4])
xline(SIC_change,'LineStyle','--')
hold off



subplot(n_panes,2,6)
%yyaxis left
area(xlabel_date,dafsd_data)
ylabel('Change in FSD cat 1 [\%]')
xline(SIC_change,'LineStyle','--')
legend(["Latg","Latm","Newi","Weld","Wave",""],'Location','northwest')
% hold off


subplot(n_panes,1,4)
%yyaxis left
area(xlabel_date,itd1_data)
ylabel('FSD SIC [\%]')
xline(SIC_change,'LineStyle','--')
legend([string(1:12),""],'Location','northwest','Orientation','horizontal','NumColumns',3)
%xlabel(xlabel_date);


