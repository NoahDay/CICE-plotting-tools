% 1. Read in JRA55 air pressure data
% 2. Read in air vectors 
% 3. For our sub-domain/region plot the air pressure and vectors
clear all 
clc
close all
addpath functions

% JRA55 filename
filename.airPressure = '/Users/noahday/GitHub/CICE-plotting-tools/observations/jra55/psl_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201901010000-201912312100.nc';
airPressure = ncread(filename.airPressure, "psl");
lat_jra = ncread(filename.airPressure, "lat");
lon_jra = ncread(filename.airPressure, "lon");

%% CICE data
clear aice air_u air_v ice_u ice_v swh  filename ice_u 
historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/2019/';
sector = "SH";
user = "noahday";
grid = 'om2';
a = dir([historydir '/*.nc']);
n_files = numel(a);
% Read in file names
for i = 1:n_files
   filename.cice(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(11:end-3);
   NFSD = ncread(filename.cice(1,:),"NFSD");
   [floe_binwidth, ~, ~, ~] = cice_parameters(NFSD);
   [aice(:,:,i), sector_mask] = data_format_sector(filename.cice(i,:), "aice",sector);
   air_u(:,:,i) = data_format_sector(filename.cice(i,:),"uatm",sector); % m/s
   air_v(:,:,i) = data_format_sector(filename.cice(i,:),"vatm",sector); % m/s
   ice_u(:,:,i) = data_format_sector(filename.cice(i,:),"uvel",sector); % m/s
   ice_v(:,:,i) = data_format_sector(filename.cice(i,:),"vvel",sector); % m/s
   %ice_u(:,:,i) = data_format(filename.cice(i,:),"uvel"); % m/s
   %ice_v(:,:,i) = data_format(filename.cice(i,:),"vvel"); % m/s
   swh(:,:,i) = data_format_sector(filename.cice(i,:),"wave_sig_ht",sector); % m
   Tair(:,:,i) = data_format_sector(filename.cice(i,:),"Tair",sector); % C
   daidtd(:,:,i) = data_format_sector(filename.cice(i,:),"daidtd",sector); % SIC
   daidtt(:,:,i) = data_format_sector(filename.cice(i,:),"daidtt",sector); % SIC
   afsdn = data_format_sector(filename.cice(i,:),"afsdn",sector); % SIC
   pancake(:,:,i) = afsdn(:,:,1,1).*floe_binwidth(1);
   dafsd_newi = data_format_sector(filename.cice(i,:),"dafsd_newi",sector); % 
   dafsd_newi_pan(:,:,i) = dafsd_newi(:,:,1).*floe_binwidth(1);
   dafsd_wave = data_format_sector(filename.cice(i,:),"dafsd_wave",sector); % 
   dafsd_wave_pan(:,:,i) = dafsd_wave(:,:,1).*floe_binwidth(1);
   dafsd_weld = data_format_sector(filename.cice(i,:),"dafsd_weld",sector); % 
   dafsd_weld_pan(:,:,i) = dafsd_weld(:,:,1).*floe_binwidth(1);
   
   ice_mask = aice(:,:,1) < 0.01;
   temp = daidtd(:,:,i);
   temp(ice_mask) = NaN;
   daidtd(:,:,i) = temp;

   temp = daidtt(:,:,i);
   temp(ice_mask) = NaN;
   daidtt(:,:,i) = temp;

   %temp = ice_u(:,:,i);
   %temp(ice_mask) = NaN;
   %ice_u(:,:,i) = temp;

   %temp = ice_v(:,:,i);
   %temp(ice_mask) = NaN;
   %ice_v(:,:,i) = temp;

   temp = pancake(:,:,i);
   temp(ice_mask) = NaN;
   pancake(:,:,i) = temp;
end
%%
i = 1;
filename.cice(i,:) = strcat(historydir,a(i).name);
dirdates(i,:) = a(i).name(11:end-3);
filedir = filename.cice(1,:);
info = ncinfo(filedir,"aice");
attributes = info.Attributes;

coord_type = "t"; % T-grid
tlat = ncread(filedir,"TLAT");
tlon = ncread(filedir,"TLON");

coord_type = "u"; % U-grid
ulat = ncread(filedir,"ULAT");
ulon = ncread(filedir,"ULON");

data_size = info.Size;
dim = length(data_size);
if data_size(1) == 320 && data_size(2) == 384
    grid = "gx1";
    row = 37;
    tlat = rearrange_matrix(tlat,row,2);
    tlon = rearrange_matrix(tlon,row,2);
    tlon = [zeros(1,384);tlon];
    tlat = [tlat(1,:); tlat];

    ulat = rearrange_matrix(ulat,row,2);
    ulon = rearrange_matrix(ulon,row,2);
    ulon = [zeros(1,384);ulon];
    ulat = [ulat(1,:); ulat];
elseif data_size(1) == 360 && data_size(2) == 300
        row = 281;
        % OM2 grid
        tlat = rearrange_matrix(tlat,row,2);
        tlon = rearrange_matrix(tlon,row,2);
        tlon = [zeros(1,300);tlon];
        tlat = [tlat(1,:); tlat];
        row = 280;
        ulat = rearrange_matrix(ulat,row,2);
        ulon = mod(rearrange_matrix(ulon,row,2)+360,360);
        ulon = [zeros(1,300);ulon];
        ulat = [ulat(1,:); ulat];
end


%% 
load('kmeans_may_sep_2019_4_classes.mat');
idx = kmeans_cluster.idx;
row_idx = kmeans_cluster.row_idx;
label_vec = kmeans_cluster.label_vec;
C = kmeans_cluster.C;
X_new = kmeans_cluster.X_new;
num_clusters = size(C,1);
region_label = {'Sheet ice','Consolidated','Wavey','Pancake'};
[len,wid] = size(tlat);

index_lat = X_new(:,end-1) == X_new(1,end-1);
index_lon = X_new(:,end) == X_new(1,end);
index_both = index_lat.*index_lon;

row_file = find(index_both);%[0,cumsum(row_idx)];
for file_number = 1:n_files-1
    row_vec = row_file(file_number)+1:row_file(file_number+1);
    file_idx = idx(row_vec);
    X_map = [file_idx, X_new(row_vec,end-1:end)];%[file_idx, X_new(:,end-1:end)];%
    k_means(:,:,file_number) = NaN.*ones(len,wid);
    for i = 1:length(file_idx)
        [lon_pos,lat_pos,~] = near2(tlon,tlat,X_map(i,3),X_map(i,2));
        k_means(lon_pos,lat_pos,file_number) = file_idx(i); 
    end
end
% w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [-90 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'on')
%     setm(w, 'parallellabel', 'off') 
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(ulat,ulon,ice_u(:,:,1))
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     colorbar
%     cmocean('balance',31)
%     caxis([0,1])
% 


% Interpolate the time data to hourly
sector = "AU";
coords = sector_coords(sector);
% September start iceh_inst.2007-07-01-75600.nc
jra_idx = 181*8 - 240 - 31*8; % May 1st (July - June - May)
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
%% AICE

sector = "SH";
coords = sector_coords(sector);
user = 'noahday';
%warning('off')
pram.ice_edge_color = 0.9*[0.4660 0.6740 0.1880];
vid_max = n_files;
line_width = 2;
FigWidth = 1500;
FigHeight = 900;
BufferWidth = 50;
Height = 100;

%pixels = 1000;
%shelf = double(isiceshelf(tlat(:,1:65),tlon(:,1:65)));
%lon_shelf = interp2(tlon(:,1:90),2);%repmat(linspace(1,365,pixels)',1,pixels);
%lat_shelf = interp2(tlat(:,1:90),2);%repmat(linspace(-90,-10,pixels),pixels,1);
%shelf = double(isiceshelf(lat_shelf,lon_shelf));
%idx = shelf == 0;
%shelf(idx) = NaN;
font_size = 18;
%shelf_vec = shelf(~idx);
%lat_shelf_vec = lat_shelf(~idx);
%lon_shelf_vec = lon_shelf(~idx);
for i = 1:vid_max
    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [100 100 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [0 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('world');
        axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        %%tightmap
        setm(w, 'Origin', [0 90 0]); 
        %setm(w, 'maplatlimit', [coords(2,1),coords(1,1)]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)]); 
        setm(w, 'maplatlimit', [-75,coords(1,1)]); setm(w, 'maplonlimit', [0,150]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',12);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat,tlon,aice(:,:,i))
        ylabel('Latitude')
        sb = scalebar('color',[1.0 1.0 1.0], 'location','sw','position',[.87 .13 .02 .85]);
        cb = colorbar; set(cb,'position',[.87 .2 .02 .65])
        cmocean('ice')
        cb.Label.String = 'Sea ice concentration'; cb.Label.Interpreter = 'latex';
        cb.FontSize = font_size+1;cb.Color = 'white';
        caxis([0,1])
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LineColor','r','LevelStep',20,'LineWidth',1.2,'ShowText','on');
        qv = quivermc(ulat,ulon,air_u(:,:,i),air_v(:,:,i),'density',50,'units','Wind speed [m/s]','reference',15,'color',[.9 .1 .1],'colormap',autumn(256),'FontSize',font_size);
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','FontSize',font_size+5)%,'Position',[1.0, -0.5, 0])
        land = shaperead('landareas', 'UseGeoCoords', true);
        %hold on
        %pcolorm(lat_shelf,lon_shelf,shelf,'FaceAlpha', 0.5)
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('aice.gif','DelayTime',0.5,'resolution',100,'overwrite',true)
    else
        gif
    end
end

%% SWH

sector = "SH";
coords = sector_coords(sector);
user = 'noahday';
pram.ice_edge_color = 0.9*[0.4660 0.6740 0.1880];
vid_max = n_files;
line_width = 2;
FigWidth = 1500;
FigHeight = 900;
BufferWidth = 50;
Height = 100;

%pixels = 1000;
%shelf = double(isiceshelf(tlat(:,1:65),tlon(:,1:65)));
%lon_shelf = interp2(tlon(:,1:90),2);%repmat(linspace(1,365,pixels)',1,pixels);
%lat_shelf = interp2(tlat(:,1:90),2);%repmat(linspace(-90,-10,pixels),pixels,1);
%shelf = double(isiceshelf(lat_shelf,lon_shelf));
%idx = shelf == 0;
%shelf(idx) = NaN;
font_size = 18;
%shelf_vec = shelf(~idx);
%lat_shelf_vec = lat_shelf(~idx);
%lon_shelf_vec = lon_shelf(~idx);
for i = 1:vid_max
    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [100 100 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [0 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('world');
        axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        %%tightmap
        setm(w, 'Origin', [0 90 0]); 
        %setm(w, 'maplatlimit', [coords(2,1),coords(1,1)]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)]); 
        setm(w, 'maplatlimit', [-75,coords(1,1)]); setm(w, 'maplonlimit', [0,150]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',12);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat,tlon,aice(:,:,i))
        ylabel('Latitude')
        sb = scalebar('color',[1.0 1.0 1.0], 'location','sw','position',[.87 .13 .02 .85]);
        cb = colorbar; set(cb,'position',[.87 .2 .02 .65])
        cmocean('haline')
        cb.Label.String = 'Significant wave height'; cb.Label.Interpreter = 'latex';
        cb.FontSize = font_size+1;cb.Color = 'white';
        caxis([0,1])
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LineColor','r','LevelStep',20,'LineWidth',1.2,'ShowText','on');
        qv = quivermc(ulat,ulon,air_u(:,:,i),air_v(:,:,i),'density',50,'units','Wind speed [m/s]','reference',15,'color',[.9 .1 .1],'colormap',autumn(256),'FontSize',font_size);
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','FontSize',font_size+5)%,'Position',[1.0, -0.5, 0])
        land = shaperead('landareas', 'UseGeoCoords', true);
        %hold on
        %pcolorm(lat_shelf,lon_shelf,shelf,'FaceAlpha', 0.5)
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('swh.gif','DelayTime',0.5,'resolution',100,'overwrite',true)
    else
        gif
    end
end


%% DAIDTD
sector = "vichi2";
coords = sector_coords(sector);
user = 'noahday';
%warning('off')
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = n_files;
line_width = 2;
FigWidth = 1920;
FigHeight = 1080;
BufferWidth = 200;
Height = 100;
%
clear aice_sector

font_size = 20;

for i = 1:vid_max
    
    [aice_sector] = sector_data(tlon,tlat,coords,aice(:,:,i));
    [tlon_sector] = sector_data(tlon,tlat,coords,tlon);
    [tlat_sector] = sector_data(tlon,tlat,coords,tlat);
    [ulon_sector] = sector_data(ulon,ulat,coords,ulon);
    [ulat_sector] = sector_data(ulon,ulat,coords,ulat);
    [daidtd_sector] = sector_data(tlon,tlat,coords,daidtd(:,:,i));
    [ice_u_sector] = sector_data(tlon,tlat,coords,ice_u(:,:,i));
    [ice_v_sector] = sector_data(tlon,tlat,coords,ice_v(:,:,i));
    %
    
    [airPressureInterp_sector] = sector_data(lon_jra,lat_jra,coords,airPressureInterp(:,:,i)/100);
    % %[lat_jra_sector] = sector_data(lon_jra,lat_jra,coords,lat_jra);
    %[lon_jra_sector] = sector_data(lon_jra,lat_jra,coords,lon_jra);
    min_lon = near1(lon_jra,coords(1,2));
    max_lon = near1(lon_jra,coords(3,2));
    
    max_lat = near1(lat_jra,coords(3,1));
    
    
    lat_jra_sector = lat_jra(1:max_lat);
    
    for j = 1:(length(lon_jra)-min_lon)+max_lon
        if j <= (length(lon_jra)-min_lon)
            lon_jra_sector(j) = lon_jra(min_lon+j-1);
        else
            lon_jra_sector(j) = lon_jra(j-(length(lon_jra)-min_lon));
        end
    
    end


    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [1 1 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [10 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('World');
        axesm eqdcylin; 
        setm(w, 'Origin', [0 28 0]); 
        setm(w, 'maplatlimit', [coords(2,1),coords(1,1)-10]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)-0]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',font_size);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat_sector,tlon_sector,daidtd_sector)
        ylabel('Latitude')
        scalebar('color',[1.0 1.0 1.0], 'location','sw');
       cb = colorbar; set(cb,'position',[.89 .15 .02 .75])
        cb.Label.String = 'Change in SIC due to dynamics [%]'; cb.Label.Interpreter = 'latex';
        cb.FontSize = font_size;cb.Color = 'white';
        cmocean('balance',31)
        caxis([-50,50])

        %[c1,h] = contourm(lat_jra_sector,lon_jra_sector,airPressureInterp_sector');%,'LevelStep',20)%,'ShowText','on')
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LevelStep',20);%,'ShowText','on')
        
        quivermc(ulat_sector,ulon_sector,ice_u_sector,ice_v_sector,'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','Position',[0.0500 -0.8663 5000],'FontSize',font_size+5)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
        

    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('daidtd.gif','DelayTime',20/n_files,'resolution',100,'overwrite',true)
    else
        gif
    end
end


%% PANCAKE
sector = "vichi2";
coords = sector_coords(sector);
user = 'noahday';
%warning('off')
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = n_files;
line_width = 1;
FigWidth = 1920;
FigHeight = 1080;
BufferWidth = 200;
Height = 100;
font_size = 20;
%
SIC = 0.15;
clear aice_sector pancake_sector ulon_sector ulat_sector lon_jra_sector



for i = 1:vid_max
    [pancake_sector] = sector_data(tlon,tlat,coords,pancake(:,:,i));
    [swh_sector] = sector_data(tlon,tlat,coords,swh(:,:,i));
    [Tair_sector] = sector_data(tlon,tlat,coords,Tair(:,:,i));
    [dafsd_newi_sector] = sector_data(tlon,tlat,coords,dafsd_weld_pan(:,:,i));
    [tlon_sector] = sector_data(tlon,tlat,coords,tlon);
    [tlat_sector] = sector_data(tlon,tlat,coords,tlat);
    [ulon_sector] = sector_data_u(ulon,ulat,coords,ulon);
    [ulat_sector] = sector_data_u(ulon,ulat,coords,ulat);
    [daidtd_sector] = sector_data(tlon,tlat,coords,daidtd(:,:,i));
    [ice_u_sector] = sector_data(tlon,tlat,coords,ice_u(:,:,i));%sector_data_u(ulon,ulat,coords,ice_u(:,:,i));
    [ice_v_sector] = sector_data(tlon,tlat,coords,ice_v(:,:,i));%sector_data_u(ulon,ulat,coords,ice_v(:,:,i));
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,tlat,tlon);
    %
    [aice_sector] = sector_data(tlon,tlat,coords,aice(:,:,i));
    icemask = aice_sector > 0.01;
    dafsd_newi_sector(~icemask) = NaN;
    [airPressureInterp_sector] = sector_data(lon_jra,lat_jra,coords,airPressureInterp(:,:,i)/100);
    min_lon = near1(lon_jra,coords(1,2));
    max_lon = near1(lon_jra,coords(3,2));
    max_lat = near1(lat_jra,coords(3,1));
    lat_jra_sector = lat_jra(1:max_lat);
    for j = 1:(length(lon_jra)-min_lon)+max_lon
        if j <= (length(lon_jra)-min_lon)
            lon_jra_sector(j) = lon_jra(min_lon+j-1);
        else
            lon_jra_sector(j) = lon_jra(j-(length(lon_jra)-min_lon));
        end
    end
    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [1 1 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [20 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('World');
        axesm eqdcylin; 
        setm(w, 'Origin', [0 28 0]); 
        setm(w, 'maplatlimit', [coords(2,1),coords(1,1)-10]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)-0]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10);  
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',font_size);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat_sector,tlon_sector,Tair_sector)
        hold on
        ylabel('Latitude')
        scalebar('color',[1.0 1.0 1.0], 'location','sw');
        cb = colorbar; set(cb,'position',[.91 .15 .02 .75])
        cb.Label.String = '10 m air temp. [C]'; cb.Label.Interpreter = 'latex';
        cb.FontSize = font_size;cb.Color = 'white';
        cmocean('thermal',31)
       
        caxis([-30,10])
        plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width) %lat_vec,lon_vec,


        %[c1,h] = contourm(lat_jra_sector,lon_jra_sector,airPressureInterp_sector');%,'LevelStep',20)%,'ShowText','on')
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LevelStep',20);%,'ShowText','on')
        
%         quivermc(ulat_sector,ulon_sector,ice_u_sector,ice_v_sector,'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
%         %quivermc(tlat,tlon,ice_u(:,:,i),ice_v(:,:,i),'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
%         %quivermc(ulat,ulon,ice_u(:,:,i),ice_v(:,:,i),'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256))'
%         t = clabelm(c1,h);
%         set(t,'Color','r'); set(t,'BackgroundColor','none'); 
%         set(t,'FontWeight','bold'); set(t,'FontSize',9);
%         title(dirdates(i,:),'Color','white','Position',[0.0500 -0.8663 5000])
%         land = shaperead('landareas', 'UseGeoCoords', true);
%         geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
%         han=axes(f,'visible','off'); 
%         han.Title.Visible='on';

        %quivermc(tlat_sector,tlon_sector,air_u_sector,air_v_sector,'density',50,'units','Wind velocity [m/s]','reference',10,'color',[.9 .1 .1],'colormap',autumn(256));
        quivermc(ulat,ulon,air_u(:,:,i),air_v(:,:,i),'density',50,'units','Wind velocity [m/s]','reference',10,'color',[.9 .1 .1],'colormap',autumn(256));
        
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','Position',[0.0500 -0.8663 5000],'FontSize',font_size+5)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
        

    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('Tair.gif','DelayTime',20/n_files,'resolution',100,'overwrite',true)
    else
        gif
    end
end

%% DAIDTT
sector = "vichi2";
coords = sector_coords(sector);
user = 'noahday';
%warning('off')
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = n_files;
line_width = 2;
FigWidth = 1920;
FigHeight = 1080;
BufferWidth = 200;
Height = 100;
%
clear aice_sector

font_size = 20;

for i = 1:vid_max
    
    [aice_sector] = sector_data(tlon,tlat,coords,aice(:,:,i));
    [tlon_sector] = sector_data(tlon,tlat,coords,tlon);
    [tlat_sector] = sector_data(tlon,tlat,coords,tlat);
    [ulon_sector] = sector_data(ulon,ulat,coords,ulon);
    [ulat_sector] = sector_data(ulon,ulat,coords,ulat);
    [daidtt_sector] = sector_data(tlon,tlat,coords,daidtt(:,:,i));
    [ice_u_sector] = sector_data(tlon,tlat,coords,ice_u(:,:,i));
    [ice_v_sector] = sector_data(tlon,tlat,coords,ice_v(:,:,i));
    %
    
    [airPressureInterp_sector] = sector_data(lon_jra,lat_jra,coords,airPressureInterp(:,:,i)/100);
    % %[lat_jra_sector] = sector_data(lon_jra,lat_jra,coords,lat_jra);
    %[lon_jra_sector] = sector_data(lon_jra,lat_jra,coords,lon_jra);
    min_lon = near1(lon_jra,coords(1,2));
    max_lon = near1(lon_jra,coords(3,2));
    
    max_lat = near1(lat_jra,coords(3,1));
    
    
    lat_jra_sector = lat_jra(1:max_lat);
    
    for j = 1:(length(lon_jra)-min_lon)+max_lon
        if j <= (length(lon_jra)-min_lon)
            lon_jra_sector(j) = lon_jra(min_lon+j-1);
        else
            lon_jra_sector(j) = lon_jra(j-(length(lon_jra)-min_lon));
        end
    
    end


    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [1 1 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [10 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('World');
        axesm eqdcylin; 
        setm(w, 'Origin', [0 28 0]); 
        setm(w, 'maplatlimit', [coords(2,1),coords(1,1)-10]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)-0]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',font_size);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat_sector,tlon_sector,daidtt_sector)
        ylabel('Latitude')
        scalebar('color',[1.0 1.0 1.0], 'location','sw');
       cb = colorbar; set(cb,'position',[.89 .15 .02 .75])
        cb.Label.String = 'Change in SIC due to thermodynamics [%]'; cb.Label.Interpreter = 'latex';
        cb.FontSize = font_size;cb.Color = 'white';
        cmocean('balance',31)
        caxis([-50,50])

        %[c1,h] = contourm(lat_jra_sector,lon_jra_sector,airPressureInterp_sector');%,'LevelStep',20)%,'ShowText','on')
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LevelStep',20);%,'ShowText','on')
        
        quivermc(ulat_sector,ulon_sector,ice_u_sector,ice_v_sector,'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','Position',[0.0500 -0.8663 5000],'FontSize',font_size+5)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
        

    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('daidtt.gif','DelayTime',20/n_files,'resolution',100,'overwrite',true)
    else
        gif
    end
end



    %% FRAME
% From Vichi "The minimum pressure found along the track simulated by ERA5 was 927"
close all
coord_grid = meshgrid(lat_jra,lon_jra);
figure
w = worldmap('world');
    %axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    %setm(w, 'Origin', [-90 0 0]);
    %setm(w, 'maplatlimit', [-90,-50]);
    %setm(w, 'maplonlimit', [0,50]);
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
    pcolorm(lat_cice,lon_cice,aice)
    c = colorbar;
    cmocean('ice')
    %pcolorm(lat_cice,lon_cice,swh)
    %c.Label.String = 'Sea level pressure [Pa]';
    caxis([0,1])
    [c1,h] = contourm(lat_jra,lon_jra,airPressure2D');%,'LevelStep',20)%,'ShowText','on')
    quivermc(lat_cice,lon_cice,air_u,air_v,'density',50,'units','Wind speed [m/s]','reference',10,'color',[.9 .1 .1])
    t = clabelm(c1,h);
    set(t,'Color','r')
    set(t,'BackgroundColor','none')
    set(t,'FontWeight','bold')
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])

    %% GIF
    % Some sample data:
    clear all
close all
t = sin(linspace(0,2*pi,30));
[X,Y,Z] = peaks(500);

% Plot the first frame:
h = surf(X,Y,Z*t(1));
shading interp
axis([-3 3 -3 3 -9 9])

% Make it fancy:
camlight
set(gca,'color','k')
set(gcf,'color','k')
caxis([min(Z(:)) max(Z(:))])

%gif('myfile.gif')
%If you want to specify certain options, include them the first time you call gif. For example, if you want a 1/24 second delay between each frame, you want the loop to run five times, and you want to use the entire figure window rather than the current axes, specifying all those options would look like this:

%gif('myfile.gif','DelayTime',1/24,'LoopCount',5)
%Or, if you want a high-resolution gif that uses export_fig, specify a resolution in units of dpi. This option is slower and creates larger files, but in some cases the difference in image quality may be significant. Here's how you might specify 400 dpi:
gif('myfile.gif','DelayTime',1/24,'resolution',400)

for k = 2:29
   set(h,'Zdata',Z*t(k))
   gif
end

web('myfile.gif')

%%

for k =1:2 
t = sin(linspace(0,2*pi,30));
[X,Y,Z] = peaks(500);

% Plot the first frame:
h = surf(X,Y,Z*t(1));
shading interp
axis([-3 3 -3 3 -9 9])

% Make it fancy:
camlight
set(gca,'color','k')
set(gcf,'color','k')
caxis([min(Z(:)) max(Z(:))])
if k==10 
    for kk = 1:48 
    gif('resolution',400) 
    end 
else 
    gif 
end
end

%%
clear all
close all
%cdt deseason
load pacific_sst
whos


sst_mean = mean(sst,3); % Along third dimension
imagescn(lon,lat,sst_mean) % n sets NAN to transparent

xlabel longitude
ylabel latitude


cmocean thermal
hold on
borders
row = near1(lat,23) % Nearest row to 23 degrees N
col = near1(lon,-115)

plot(lon(col),lat(row),'ko')

sst1 = squeeze(sst(row,col,:));
%%
clf

plot(t,sst1)
ylabel temperate
xlabel time
datetick

sst_ds = deseason(sst1,t);
hold on
plot(t,sst_ds)
polyplot(t,sst1) % Linear trend, least squares fit
polyplot(t,sst_ds)
trend(sst1,12) % Trend per year, 12 months
%%
clf 
sst_tr = trend(sst,12);
imagescn(lon,lat,sst_tr)
cb = colorbar;
ylabel(cb,'sst trend \circC/ye')
cmocean('balance','pivot')

[sst_tr,p]  = trend(sst,12); % get p values
mk = mann_kendall(sst); % is this statistically significant 1 = yes
hold 
stipple(lon,lat,mk)


%% Calculating curl


clear aice air_u air_v ice_u ice_v swh  filename ice_u dirdates
historydir = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history_h/';
sector = "SH";
user = "noahday";
grid = 'om2';

a = dir([historydir '/*.nc']);
n_files = numel(a);

% Read in file names
for i = 1:n_files
    filename.cice(i,:) = strcat(historydir,a(i).name);
    dirdates(i,:) = a(i).name(11:end-3);
end
i = 60;
NFSD = ncread(filename.cice(1,:),"NFSD");
[lat,lon,~,ulat,ulon] = grid_read('om2');

[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
[aice(:,:,i), sector_mask] = data_format_sector(filename.cice(i,:), "aice",sector);
air_u(:,:,i) = data_format_sector(filename.cice(i,:),"uatm",sector); % m/s
air_v(:,:,i) = data_format_sector(filename.cice(i,:),"vatm",sector); % m/s
ice_u(:,:,i) = data_format_sector(filename.cice(i,:),"uvel",sector); % m/s
ice_v(:,:,i) = data_format_sector(filename.cice(i,:),"vvel",sector); % m/s

%cav = curl(lon,lat,ice_u,ice_v);

% pcolor(lon,lat,cav); shading interp
% ylim([-90,-60])
% colorbar
% hold on; quiver(lon,lat,ice_u,ice_v)

[Lat,Lon] = meshgrid(lat(1,:),lon(:,180));
[Lat,Lon,u,v] = recenter(Lat,Lon,ice_u(:,:,i),ice_v(:,:,i));
C = cdtcurl(Lat,Lon,u,v);
C = C.*sign(Lat); % Multiply in the south with -1
%C = cdtcurl(lat',lon',ice_u',ice_v');
%
lat_pos = near1(Lat(1,:),-62);
lon_pos = near1(Lon(:,180), 8);

C(lon_pos,lat_pos) = 10;

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
    pcolorm(Lat,Lon,C)
    ylabel('Latitude')
    scalebar('color',[0.0 0.0 0.0], 'location','sw');
    cb = colorbar; set(cb,'position',[.89 .15 .02 .75])
    cb.Label.String = 'CURL '; cb.Label.Interpreter = 'latex';
    cb.FontSize = font_size;cb.Color = 'black';
    cmocean('balance',31)
    caxis([-10^(-5),10^(-5)])
    [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LevelStep',20);%,'ShowText','on')
    quivermc(ulat,ulon,ice_u(:,:,i),ice_v(:,:,i),'density',50,'units','Ice drift [m/s]','reference',0.2,'color',[.9 .1 .1],'colormap',cool(256));
    t = clabelm(c1,h);
    set(t,'Color','r'); set(t,'BackgroundColor','none'); 
    set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
    %title(dirdates(i,:),'Color','black','Position',[0.0500 -0.8663 5000],'FontSize',font_size+5)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])



%%
ice_mag = sqrt(ice_u(:,:,i).^2 + ice_v(:,:,i).^2);
ice_mag = reshape(ice_mag,[numel(ice_mag),1]);
C_vec = reshape(C,[numel(C),1]);
figure
scatter(C_vec,ice_mag)
%%
close all
 load wind
   k = 4; 
   x = x(:,:,k); % LON
   y = y(:,:,k); % LAT
   u = u(:,:,k); 
   v = v(:,:,k); 
 Cz = cdtcurl(y,x,u,v)  
   cav = curl(x,y,u,v);
   pcolor(x,y,cav); shading interp
   hold on; quiver(x,y,u,v)

%% k-means with cyclone data (k_means)

sector = "SH";
coords = sector_coords(sector);
user = 'noahday';
pram.ice_edge_color = 0.9*[0.4660 0.6740 0.1880];
vid_max = n_files-1;
line_width = 2;
FigWidth = 1500;
FigHeight = 900;
BufferWidth = 50;
Height = 100;
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];

font_size = 18;
SIC = 0.15;
region_label = {'Sheet ice','Consolidated','Wavey','Pancake'};
for i = 1:vid_max
    [tlon_sector] = sector_data(tlon,tlat,coords,tlon);
    [tlat_sector] = sector_data(tlon,tlat,coords,tlat);
    [ulon_sector] = sector_data_u(ulon,ulat,coords,ulon);
    [ulat_sector] = sector_data_u(ulon,ulat,coords,ulat);
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,tlat,tlon);

    f = figure('PaperPositionMode','manual','visible','off');
    f.Position = [100 100 540 200];
    AX = gca;
    f.Resize = 'off'; f.Units = 'points'; f.PaperUnits = 'points';
    AX.Units = 'points'; AX.Clipping = 'on'; AX.PositionConstraint = 'innerposition';
    AX.InnerPosition = [0 0 FigWidth-BufferWidth FigHeight-Height]*72/96; % converting from pixels to points
    f.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
    f.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
    camlight
    set(gca,'color','k')
    set(gcf,'color','k')
    w = worldmap('world');
        axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        %%tightmap
        setm(w, 'Origin', [0 90 0]); 
        %setm(w, 'maplatlimit', [coords(2,1),coords(1,1)]); setm(w, 'maplonlimit', [coords(1,2)-360,coords(3,2)]); 
        setm(w, 'maplatlimit', [-75,coords(1,1)]); setm(w, 'maplonlimit', [0,150]); 
        setm(w, 'meridianlabel', 'on'); setm(w, 'parallellabel', 'on'); 
        setm(w, 'mlabellocation', 30); setm(w, 'plabellocation', 10); 
        setm(w, 'mlabelparallel', -80,'FontColor','white','FontSize',12);
        setm(w, 'grid', 'on'); setm(w, 'labelrotation', 'on');
        pcolorm(tlat,tlon,k_means(:,:,i))
        ylabel('Latitude')
        sb = scalebar('color',[1.0 1.0 1.0], 'location','sw','position',[.87 .13 .02 .85]);
        cb = colorbar; set(cb,'position',[.87 .2 .02 .65])
        cmocean('deep',num_clusters)
        %cb.Label.String = 'Sea ice concentration'; 
        cb.FontSize = font_size+1;cb.Color = 'white';
        cb.TickLabels = region_label;
        cb.Ticks = 1:num_clusters;
        cb.Label.Interpreter = 'latex';
        caxis([1,num_clusters])
        [c1,h] = contourm(lat_jra,lon_jra,airPressureInterp(:,:,i)'/100,'LineColor','r','LevelStep',20,'LineWidth',1.2,'ShowText','on');
        qv = quivermc(ulat,ulon,air_u(:,:,i),air_v(:,:,i),'density',50,'units','Wind speed [m/s]','reference',20,'color',[.9 .1 .1],'colormap',autumn(256),'FontSize',font_size);
        t = clabelm(c1,h);
        set(t,'Color','r'); set(t,'BackgroundColor','none'); 
        set(t,'FontWeight','bold'); set(t,'FontSize',font_size);
        title(dirdates(i,:),'Color','white','FontSize',font_size+5)%,'Position',[1.0, -0.5, 0])
        land = shaperead('landareas', 'UseGeoCoords', true);
        plotm(lat_ice_edge,lon_ice_edge+1,'-','color',pram.ice_edge_color,'LineWidth',line_width)
        geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
        han=axes(f,'visible','off'); 
        han.Title.Visible='on';
        

    set(gcf,'Color',[0 0 0]); % color of the frame around the figure
    set(gca,'Color','k')%color for the plot area
    set(gca,'XColor',[1 1 1]); % Set RGB value to what you want
    set(gca,'YColor',[1 1 1]); % Set RGB value to what you want
    if i == 1
        gif('kmeans_cyclone.gif','DelayTime',0.5,'resolution',100,'overwrite',true)
    else
        gif
    end
end





%% Functions

function [aice_sector] = sector_data(tlon,tlat,coords,aice)
    [len,wid] = size(tlon);
    min_lat = near1(tlon(:,1),coords(1,2)); % lat should be lon
    max_lat = near1(tlon(:,1),coords(3,2));
    
    if size(tlat(1,:)) == 1
        max_height = near1(tlat,coords(3,1));
    else
        max_height = near1(tlat(1,:),coords(3,1));
    end
    if min_lat > max_lat
        % Take the edge slices
        for i = 1:len-min_lat + max_lat
            if i <= len-min_lat
                % West slice
                aice_sector(i,:) = aice(min_lat+i-1,1:max_height,1);
            else
                % East slice
                aice_sector(i,:) = aice(i-(len-min_lat),1:max_height,1);
            end
        end
    
    else
        % Take middle bit
        for i = 1:max_lat-min_lat
            %disp(size(aice(min_lat+i-1,1:max_height)))
            aice_sector(i,:) = aice(min_lat+i-1,1:max_height);
        end
    end

end


function [aice_sector] = sector_data_u(tlon,tlat,coords,aice)
    [len,wid] = size(tlon);
    min_lat = 346;%near1(tlon(:,1),coords(1,2));
    max_lat = 41;%near1(tlon(:,1),coords(3,2));
    
    if size(tlat(1,:)) == 1
        max_height = near1(tlat,coords(3,1));
    else
        max_height = near1(tlat(1,:),coords(3,1));
    end
    
    if min_lat > max_lat
        % Take the edge slices
        for i = 1:len-min_lat + max_lat
            if i <= len-min_lat
                % West slice
                aice_sector(i,:) = aice(min_lat+i-1,1:max_height,1);
            else
                % East slice
                aice_sector(i,:) = aice(i-(len-min_lat),1:max_height,1);
            end
        end
    
    else
        % Take middle bit
        
        for i = 1:max_lat-min_lat
            aice_sector(i) = aice(min_lat+i-1,1:max_height);
        end
    end

end
