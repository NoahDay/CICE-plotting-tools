%% Case info
clear all
clc
close all
addpath functions

historydir = '/Users/noahday/Maths1/spectrum_fixed/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end

% Get the data
SIC = 0.15;

[lat,lon] = grid_read("gx1");
NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
sector = "SA";
for i = 1:n_files
    [dafsd.latm, sector_mask] = data_format_sector(filenames(i,:),"dafsd_latm",sector);
    dafsd.latg = data_format_sector(filenames(i,:),"dafsd_latg",sector);
    dafsd.weld = data_format_sector(filenames(i,:),"dafsd_weld",sector);
    dafsd.newi = data_format_sector(filenames(i,:),"dafsd_newi",sector);
    dafsd.wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
    aice(:,:,i)  = data_format_sector(filenames(i,:),"aice",sector);
    afsd_data = data_format_sector(filenames(i,:),"afsd",sector);
    swh(:,:,i) = data_format_sector(filenames(i,:),"wave_sig_ht",sector);
    idx = aice(:,:,i) > eps; 
    for nf = 1:length(floe_binwidth)
        temp(:,:) = dafsd.latg(:,:,nf);
        afsd.latg(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.latm(:,:,nf);
        afsd.latm(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.newi(:,:,nf);
        afsd.newi(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.weld(:,:,nf);
        afsd.weld(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.wave(:,:,nf);
        afsd.wave(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
    end

end

[len,wid, ~] = size(aice);
for k = 1:n_files
    for i = 1:len
        for j = 1:wid
            afsd.latg2d(i,j,k) = sum(afsd.latg(i,j,:,k));
            afsd.latm2d(i,j,k) = sum(afsd.latm(i,j,:,k));
            afsd.newi2d(i,j,k) = sum(afsd.newi(i,j,:,k));
            afsd.weld2d(i,j,k) = sum(afsd.weld(i,j,:,k));
            afsd.wave2d(i,j,k) = sum(afsd.wave(i,j,:,k));
        end
    end
end
%% Plotting
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
line_width = 1;
close all
clear writerObj
video_name = strcat("spec_fixed", '_', "dafsd", '_', '2022_08_15', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
set(writerObj,'FrameRate',1); % 0.5 = 2 seconds per frame
%writerObj.Quality = 70;
user = 'noahday';
% open the writer
open(writerObj);
% Plotting and saving
clear C2
month_strings = ["Jan. 2005","Feb. 2005","Mar. 2005","Apr. 2005","May 2005","June 2005","July 2005","Aug. 2005","Sep. 2005","Oct. 2005","Nov. 2005","Dec. 2005","Jan. 2006","Feb. 2006","Mar. 2006","Apr. 2006","May 2006","June 2006","July 2006","Aug. 2006","Sep. 2006","Oct. 2006","Nov. 2006","Dec. 2006","Jan. 2007","Feb. 2007","Mar. 2007","Apr. 2007","May 2007","June 2007","July 2007","Aug. 2007","Sep. 2007","Oct. 2007","Nov. 2007","Dec. 2007","Jan. 2008","Feb. 2008","Mar. 2008","Apr. 2008","May 2008","June 2008","July 2008","Aug. 2008","Sep. 2008","Oct. 2008","Nov. 2008","Dec. 2008","Jan. 2009","Feb. 2009","Mar. 2009","Apr. 2009","May 2009","June 2009","July 2009","Aug. 2009","Sep. 2009","Oct. 2009","Nov. 2009","Dec. 2009"];
vid_max = 36;
for i = 24:vid_max
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice(:,:,i),SIC,sector,lat,lon);
    basicwaitbar(i,12,"Plotting and saving")
    
   % Plot the map
   conFigure(30,1.1)
   %set(gcf, 'Position',  [0, 0, 1782, 1830])
   f = figure;
       subplot(1,7,1)
       cmocean('ice',11)
       [p,a] = map_plot(aice(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       cmocean('ice',11)
       caxis([0,1])
       a.Label.String = " ";
       a.Label.Interpreter = 'latex';
       title('Aice','Interpreter','latex')

       subplot(1,7,2)
       [p,a] = map_plot(afsd.latg2d(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       C = cmocean('balance',11);
       C2 = C;
       colormap(C2)
       caxis([-10^(-3),10^(-3)])
       a.Label.String = " ";
       a.Label.Interpreter = 'latex';
       title('Lat. growth','Interpreter','latex')
       

       subplot(1,7,3)
       [p,a] = map_plot(afsd.latm2d(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       C = cmocean('balance',11);
       C2 = C;
       colormap(C2)
       caxis([-10^(-1),10^(-1)])
       a.Label.String = " ";
       a.Label.Interpreter = 'latex';
       title('Lat. melt','Interpreter','latex')

       subplot(1,7,4)
       [p,a] = map_plot(afsd.newi2d(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       C = cmocean('balance',11);
       C2 = C;
       colormap(C2)
       caxis([-10^(-1),10^(-1)])
       a.Label.String = " ";
       a.Label.Interpreter = 'latex';
       title('New ice','Interpreter','latex')

       subplot(1,7,5)
       [p,a] = map_plot(afsd.weld2d(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       C = cmocean('balance',11);
       C2 = C;
       colormap(C2)
       caxis([-10^(-1),10^(-1)])
       a.Label.String = " ";
       a.Label.Interpreter = 'latex';
       title('Welding','Interpreter','latex')

       subplot(1,7,6)
       [p,a] = map_plot(afsd.wave2d(:,:,i),"aice",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       C = cmocean('balance',11);
       C2 = C;
       colormap(C2)
       caxis([-10^(-1),10^(-1)])
       a.Label.String = "Change in AFSD [$\%$]";
       a.Label.Interpreter = 'latex';
       title('Wave','Interpreter','latex')

       subplot(1,7,7)
       [p,a] = map_plot(swh(:,:,i),"wave_sig_ht",sector);
       plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
       %C = cmocean('dense',11);
       %C2 = C;
       colormap(C2)
       caxis([-3,3])
       a.Label.String = "Sig. wave height [m]";
       a.Label.Interpreter = 'latex';
       title('SWH','Interpreter','latex')


    han=axes(f,'visible','off'); 
    han.Title.Visible='on';
    %han.XLabel.Visible='on';
    %han.YLabel.Visible='on';
    %ylabel(han,'Change in SIC [\%/day]');
    %xlabel(han,'Floe size categories');
    title(han,month_strings(i));
    f.Position = [100 100 540 200];
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %exportgraphics(f,figname,'ContentType','vector')
    saveas(f,fullfile(filedir, figname));
    
    close(f)
    close all
end    
       
% iterate over each image
 for k=24:vid_max
     % use imread to read the image
     figname = sprintf('frames/image%d.png',k);
     img = imread(figname);
     % resize the imagei_vec
     %img = imresize(img,2);
     % convert the image to a frame using im2frame
     frame = im2frame(img);
     % write the frame to the video
     %[height width ~] = size(frame.cdata);
     %writerObj.Height = height;
     %writerObj.Width = width;
     writeVideo(writerObj,frame);
 end

 % close the writer
 close(writerObj);
