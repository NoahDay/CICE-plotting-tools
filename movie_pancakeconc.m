% movie_pancakeconc is a script that will produce movies of pancake ice
% given a defintion.

% Definition 1: Areal ice concentration in FSD cat 1 and ITD cat 1.
% Definition 2: Areal ice concentration in FSD cat 1 and ITD cat 1/ aice.
%% Case info
clear all
clc
close all
addpath functions

pancake_def = [2];
plot_type = "univariate"; % Just plot a single def
% bivariate takes the difference pancake_def(1)-pancake_def(2);
historydir = '/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/';
sector = "SH";
user = "noahday";
%'/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/';
% '/Volumes/NoahDay5TB/cases/monthwim/history/';
a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end


% Get the data
SIC = 0.01;
pancakeconc(:,:,:)= get_pancakes(pancake_def,n_files,filenames,SIC,plot_type);


% Plotting
clear writerObj
video_name = strcat("pancake_def2", '_', "monthwim", '_', '2022_06_15', '.mp4');
writerObj = VideoWriter(video_name,'MPEG-4');
set(writerObj,'FrameRate',n_files/20); % 0.5 = 2 seconds per frame
%writerObj.Quality = 70;

% open the writer
open(writerObj);
% Plotting and saving
clear C2
month_strings = ["Jan. 2005","Feb. 2005","Mar. 2005","Apr. 2005","May 2005","June 2005","July 2005","Aug. 2005","Sep. 2005","Oct. 2005","Nov. 2005","Dec. 2005","Jan. 2006","Feb. 2006","Mar. 2006","Apr. 2006","May 2006","June 2006","July 2006","Aug. 2006","Sep. 2006","Oct. 2006","Nov. 2006","Dec. 2006","Jan. 2007","Feb. 2007","Mar. 2007","Apr. 2007","May 2007","June 2007","July 2007","Aug. 2007","Sep. 2007","Oct. 2007","Nov. 2007","Dec. 2007","Jan. 2008","Feb. 2008","Mar. 2008","Apr. 2008","May 2008","June 2008","July 2008","Aug. 2008","Sep. 2008","Oct. 2008","Nov. 2008","Dec. 2008","Jan. 2009","Feb. 2009","Mar. 2009","Apr. 2009","May 2009","June 2009","July 2009","Aug. 2009","Sep. 2009","Oct. 2009","Nov. 2009","Dec. 2009"];
for i = 1:n_files
    close all
   % Plot the map
   conFigure(30,1.1)
   %set(gcf, 'Position',  [0, 0, 1782, 1830])
   f = figure;


   if plot_type == "univariate"
       [p,a] = map_plot(pancakeconc(:,:,i),"aice",sector);
       %C = colormap(cool(20));
       C = cmocean('thermal',10);
       C2 = C;
       % For 'cool'
       %C2(:,1) = sqrt(C(:,1));
       %C2(:,2) = (C(:,2)).^2;
       %C2(:,3) = sqrt(C(:,3));
       %C2(1,:) = ones(1,3); % 0.05 and less is considered white!
       %C2(2,:) = C(1,:);
       %C2(2:end,:) = C(1:end-1,:);
       colormap(C2)
   else
       [p,a] = map_plot(pancakeconc(:,:,i),"aice",sector,"gx1",[-0.5,0.5]);
       C = cmocean('balance',20);
       colormap(C)
   end
   a.Label.String = "Pancake ice concentration";
   title(month_strings(i),'Interpreter','latex')
    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %exportgraphics(f,figname,'ContentType','vector')
    saveas(f,fullfile(filedir, figname));
end    
       
% iterate over each image
 for k=1:n_files
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


 %% Ice thickness

aice = data_format_sector(filenames(21,:),"aice",sector);
ice_mask = aice > SIC;
hi(:,:) = data_format_sector(filenames(21,:),"hi",sector);
% Define pancakes as the thinnest NCAT and smallest NFSD
temp(:,:) =  hi;
% apply mask
temp(~ice_mask) = NaN;
hi(:,:) = temp;
% Plot the map
conFigure(11,1.1)
%set(gcf, 'Position',  [0, 0, 1782, 1830])
f = figure;
[p,a] = map_plot(hi(:,:),"hi",sector,"gx1",[0,3]);
a.Label.String = "Ice thickness [m]";
exportgraphics(f,'hi.pdf','ContentType','vector')


%% FSDrad

[lat lon] = grid_read("gx1");
aice = data_format_sector(filenames(21,:),"aice",sector);
ice_mask = aice > SIC;
hi(:,:) = data_format_sector(filenames(21,:),"fsdrad",sector);
% Define pancakes as the thinnest NCAT and smallest NFSD
temp(:,:) =  hi;
% apply mask
temp(~ice_mask) = NaN;
fsdraddata(:,:) = temp;
% Plot the map
% conFigure(11,1.1)
% %set(gcf, 'Position',  [0, 0, 1782, 1830])
% f = figure;
% [p,a] = map_plot(hi(:,:),"fsdrad",sector,"gx1",[0,3000]);
% colors = cmocean('thermal');
% colormap(colors(end:-1:1,:))
% a.Label.String = "Mean floe size [m]";
% exportgraphics(f,'ra.pdf','ContentType','vector')


close all
conFigure(11,1.1)
f = figure;
%[p,a] = map_plot(fsdraddata,"fsdrad",sector);
w = worldmap('world');
            axesm eqaazim; %, wetch
            setm(w, 'Origin', [-90 0 0]);
            setm(w, 'maplatlimit', [-90,-55]);
            setm(w, 'maplonlimit', [-180,-55]);
            setm(w, 'meridianlabel', 'on')
            setm(w, 'parallellabel', 'off')
            setm(w, 'mlabellocation', 60);
            setm(w, 'plabellocation', 10);
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlinelimit', [-75 -55]);
            setm(w, 'plinelimit', [-75 -55]);
            setm(w, 'grid', 'off');
            %setm(w, 'frame', 'on');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,fsdraddata)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
            a = colorbar;
            s = scaleruler;
            setm(handlem('scaleruler'), ...
            'XLoc',0.3, ... '
            'YLoc',-0.52, ...
            'TickDir','down', ...
            'MajorTick',0:1000:1000, ...
            'MinorTick',0:500:500, ...
            'MajorTickLength',km2nm(150),...
            'MinorTickLength',km2nm(150))

a.TickLabelInterpreter = 'latex';
a.Label.Interpreter = 'latex';
a.Label.String = 'Mean floe size radius (m)';

w.ZTickLabel = 'test';
w.ColorScale = 'log';
w.FontName = 'CMU Serif';
cmocean('matter');
%colormap parula
exportgraphics(f,'fsdradsep.pdf','ContentType','vector')

%% Functions
function [data] = get_pancakes(definition,n_files,filenames,SIC,plot_type)
    sector = "SH"; % Southern Hemisphere
    if plot_type == "univariate"
        if definition == 1
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1);
                % apply mask
                temp(~ice_mask) = NaN;
                data(:,:,i) = temp;
            end
        elseif definition == 2
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1)./aice(:,:);
                % apply mask
                temp(~ice_mask) = NaN;
                data(:,:,i) = temp;
            end
        else 
            error('No pancake definition specified')
        end
    elseif plot_type == "bivariate"
        % Calcualte def 1
        if definition(1) == 1
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1);
                % apply mask
                temp(~ice_mask) = NaN;
                temp1(:,:,i) = temp;
            end
        elseif definition(1) == 2
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1)./aice(:,:);
                % apply mask
                temp(~ice_mask) = NaN;
                temp1(:,:,i) = temp;
            end
            
        end
        if definition(2) == 1 % Calculate def 2
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1);
                % apply mask
                temp(~ice_mask) = NaN;
                temp2(:,:,i) = temp;
            end
        elseif definition(2) == 2
            for i = 1:n_files
                NFSD = ncread(filenames(1,:),"NFSD");
                aice = data_format_sector(filenames(i,:),"aice",sector);
                ice_mask = aice > SIC;
                afsdn(:,:,:,:) = data_format_sector(filenames(i,:),"afsdn",sector);
                % Define pancakes as the thinnest NCAT and smallest NFSD
                temp(:,:) =  afsdn(:,:,1,1).*NFSD(1)./aice(:,:);
                % apply mask
                temp(~ice_mask) = NaN;
                temp2(:,:,i) = temp;
            end
        end
        data(:,:,:) = temp1 - temp2;
    else
        error('No plot_type specified.')
    end
     
end