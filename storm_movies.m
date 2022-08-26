clear all
clc
close all
addpath functions

historydir = '/Users/noahday/Maths1/spectrum_fixed_h/history/';
sector = "vichi";
user = "noahday";
%'/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/';
% '/Volumes/NoahDay5TB/cases/monthwim/history/';
a = dir([historydir '/*.nc']);
n_files = numel(a);

% Read in file names
for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end
NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");

filedir = filenames(1,:);
variable = "aice";
% Get the info
info = ncinfo(filedir,variable);
attributes = info.Attributes;
coord_att = attributes(3); % Extract coordinate info
coord_string = coord_att.Value;
if coord_string(1:9) == 'TLON TLAT'
    coord_type = "t"; % T-grid
    coord_type = "t"; % T-grid
    lat = ncread(filedir,"TLAT");
    lon = ncread(filedir,"TLON");
elseif coord_string(1:9) == 'ULON ULAT'
    coord_type = "u"; % U-grid
    lat = ncread(filedir,"ULAT");
    lon = ncread(filedir,"ULON");
else
    error('ND: Unspecified grid type');
end
data_size = info.Size;
idx = (data_size == 1);
dim = length(data_size(~idx));
if data_size(1) == 320 && data_size(2) == 384
    grid = "gx1";
    row = 37;
    lat = rearrange_matrix(lat,row,2);
    lon = rearrange_matrix(lon,row,2);

    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
else
    [lat,lon,row] = grid_read('om2');
    %error('ND: Unrecognised grid size.')
end
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);

%%

% Get this data
SIC = 0.15;
sector = "SH";
for i = 1:n_files
    temp = data_format_sector(filenames(i,:),"aice",sector);
    ice_mask(:,:,i) = temp > 0.01;
    %temp_mask = temp > 0.01;
    temp(~ice_mask(:,:,i)) = NaN;
    aice(:,:,i) = temp;
    air_temp(:,:,i) = data_format_sector(filenames(i,:),"Tair",sector); % Celsius
    air_u(:,:,i) = data_format_sector(filenames(i,:),"uatm",sector); % m/s
    air_v(:,:,i) = data_format_sector(filenames(i,:),"vatm",sector); % m/s
    ice_u(:,:,i) = data_format_sector(filenames(i,:),"uvel",sector); % m/s
    ice_v(:,:,i) = data_format_sector(filenames(i,:),"vvel",sector); % m/s
    daidtd(:,:,i) = data_format_sector(filenames(i,:),"daidtd",sector); % area tendency dynamics, %/day
    daidtt(:,:,i) = data_format_sector(filenames(i,:),"daidtt",sector); % area tendency thermo,%/day
    dvidtt(:,:,i) = data_format_sector(filenames(i,:),"dvidtt",sector); % volume tendency thermo, %/day
    dvidtd(:,:,i) = data_format_sector(filenames(i,:),"dvidtd",sector); % volume tendency dynamimcs, %/day
    sst(:,:,i) = data_format_sector(filenames(i,:),"sst",sector); % sea suface temp, degrees C
    fsdrad(:,:,i) = data_format_sector(filenames(i,:),"fsdrad",sector); % mean floe size m
    wave_sig_ht(:,:,i) = data_format_sector(filenames(i,:),"wave_sig_ht",sector); % SWH m
    afsdn(:,:,:,:,i) = data_format_sector(filenames(i,:),"afsdn",sector); % SWH m
    pancake(:,:,i) = afsdn(:,:,1,1,i).*floe_binwidth(1);
end


%% map

sector = "SH";
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = 28;
line_width = 2;
i = 1;
conFigure(15,3)
f = figure('PaperPosition',[.25 .25 50 30]);
[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
conFigure(15,3)
[w, a, output_data] = map_plot(wave_sig_ht(:,:,i),"aice",sector); % 7 Means average from July

quivermc(lat,lon,ice_u(:,:,i),ice_v(:,:,i),'density',100,'units','Ice drift [m/s]','reference',0.2)
plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width) %lat_vec,lon_vec,
cmocean('amp',10)
a.Label.String = "$H_s$ [$m$]";
a.Label.Interpreter = "latex";
caxis([0,3]);




%% SWH

sector = "vichi";
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = 28;
line_width = 2;
for i = 1:vid_max
    conFigure(15,3)
    f = figure('PaperPosition',[.25 .25 50 30]);
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
    conFigure(15,3)
    [w, a, output_data] = map_plot(wave_sig_ht(:,:,i),"aice","vichi"); % 7 Means average from July
    
    quivermc(lat,lon,ice_u(:,:,i),ice_v(:,:,i),'density',100,'units','Ice drift [m/s]','reference',0.2)
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width) %lat_vec,lon_vec,
    cmocean('haline',10)
    a.Label.String = "$H_s$ [$m$]";
    a.Label.Interpreter = "latex";
    caxis([0,3]);
    
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %exportgraphics(f,figname,'ContentType','vector')
    saveas(f,fullfile(filedir, figname));
    close all
end

clear writerObj
video_name = strcat("spec_fixed_h", '_', "swh", '_', '2022_08_22', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
set(writerObj,'FrameRate',1); % 0.5 = 2 seconds per frame
%writerObj.Quality = 70;
user = 'noahday';
% open the writer
open(writerObj);
% iterate over each image
 for k=1:vid_max
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     frame = im2frame(img);
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);


 %% aice

sector = "vichi";
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
vid_max = 28;
line_width = 2;
for i = 1:vid_max
    conFigure(15,3)
    f = figure('PaperPosition',[.25 .25 50 30]);
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
    conFigure(15,3)
    [w, a, output_data] = map_plot(aice(:,:,i),"aice","vichi"); % 7 Means average from July
    
    quivermc(lat,lon,ice_u(:,:,i),ice_v(:,:,i),'density',100,'units','Ice drift [m/s]','reference',0.2)
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width) %lat_vec,lon_vec,
    cmocean('ice',10)
    a.Label.String = "$H_s$ [$m$]";
    a.Label.Interpreter = "latex";
    caxis([0,1]);
    
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %exportgraphics(f,figname,'ContentType','vector')
    saveas(f,fullfile(filedir, figname));
    close all
end

clear writerObj
video_name = strcat("spec_fixed_h", '_', "swh", '_', '2022_08_22', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
set(writerObj,'FrameRate',1); % 0.5 = 2 seconds per frame
%writerObj.Quality = 70;
user = 'noahday';
% open the writer
open(writerObj);
% iterate over each image
 for k=1:vid_max
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     frame = im2frame(img);
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);

    %% pancake ice vel
sector = "vichi";
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 2;
vid_max = 28;
for i = 1:vid_max
     conFigure(15,3)
    f = figure('PaperPosition',[.25 .25 50 30]);
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
    [w, a, output_data] = map_plot(pancake(:,:,i)./aice(:,:,i),"fsdrad","vichi"); % 7 Means average from July
    
    quivermc(lat,lon,ice_u(:,:,i),ice_v(:,:,i),'density',100,'units','Ice speed [m/s]','reference',0.2)
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width) %lat_vec,lon_vec,
    cmocean('thermal',11)
    a.Label.String = "FSTD(1,1)/aice [%]";
    a.Label.Interpreter = "latex";
    caxis([0,1]);
    
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %exportgraphics(f,figname,'ContentType','vector')
    saveas(f,fullfile(filedir, figname));
     close all
end



clear writerObj
video_name = strcat("spec_fixed_h", '_', "pancake", '_', '2022_08_22', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
set(writerObj,'FrameRate',1); % 0.5 = 2 seconds per frame
%writerObj.Quality = 70;
user = 'noahday';
% open the writer
open(writerObj);
% iterate over each image
 for k=1:vid_max
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     frame = im2frame(img);
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);
