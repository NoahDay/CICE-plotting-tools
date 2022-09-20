%Initial_ice is a report of the ice initial conditions for CICE. 
%   It is taken from point observations in Perovich et al. (2014) from the
%   Arctic. This investigation will be targetted at the FSD, ITD, ice 
%   extent around Antarctica.
%% Read in the data.
clear
close all
addpath functions
filedir = [%"cases/init/history/iced_gx1_v6.2005-01-01.nc";
           "cases/init/history/iceh_ic.2005-01-01-03600.nc";
           %"cases/init/history/iceh_inst.2005-01-01-03600.nc";
           "cases/init/history/iceh.2005-01-01.nc";
           "cases/init/history/iceh.2005-09-30.nc";
           "cases/init/history/iceh.2005-12-31.nc"];
           
filecase = [%"iced\_gx1\_v6.2005-01-01";
            "iceh\_ic.2005-01-01-03600";
            %"iceh\_inst.2005-01-01-03600";
            "iceh.2005-01-01";
            "iceh.2005-09-30";
            "iceh.2005-12-31.nc"];
% Files
% Initial conditions forcing: iced_gx1_v6.2005-01-01.nc
% Initial conditions of CICE: iceh_ic.2005-01-01-03600.nc
% First timestep (1 hour): iceh_inst.2005-01-01-03600.nc
% After one day: iceh.2005-01-01.nc
ncdisp(filedir(1))

% Parameters
fstd_switch = 0;
coords = [-45,20;-65,20;-45,40;-65,40]; %(NW;SW,NE,SE)
grid = 'gx1';
[lat,lon,row] = grid_read(grid);
%% iced_gx1_v6.2005-01-01.nc

data = data_format(filedir(1),'aicen',row,lat,lon,3); %aicen
t=tiledlayout(2,2);
figure(1)
nexttile
norm_ITD = frac_pdf(data,0);
x_axis = 1:5;%linspace(0,5,length(final_data));
area(x_axis,norm_ITD)
title("Initial ITD provided from iced\_gx1\_v6.2005-01-01")
ylabel("probability")
xlabel("ITD bins")

data = data(:,:,1);
color_map = seaicecolormap();
latitude = [-90,-60];
longitude = [10,50];
nexttile
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-40]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabelparallel', 0);
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    title("AICE category 1",'interpreter','latex','FontSize', 12)
    
% SST
nexttile
data = data_format(filedir(1),'sst',row,lat,lon,3); %aicen
w = worldmap('world');
    axesm eqaazim; 
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-40]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabelparallel', 0);
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    title("SST",'interpreter','latex','FontSize', 12)
    
    
%% iceh_ic.2005-01-01-03600.nc    
for j = 2:length(filedir)
    figure(j)
    t=tiledlayout(2,2);
    nexttile
    ncdisp(filedir(j))    
    data = data_format(filedir(j),'aicen',row,lat,lon,3); %aicen
    % ITD
    norm_ITD = frac_pdf(data);
    x_axis = 1:5;%linspace(0,5,length(final_data));
    area(x_axis,norm_ITD)
    title(strcat("ITD from ", filecase(j)))
    ylabel("probability")
    xlabel("ITD bins")
    nexttile
    % FSD
    data = data_format(filedir(j),'afsd_d',row,lat,lon,3); %aicen
    norm_FSD = frac_pdf(data);
    x_axis = 1:length(norm_FSD);
    area(x_axis,norm_FSD)
    title(strcat("FSD from ", filecase(j)))
    ylabel("probability")
    xlabel("FSD bins")
    nexttile
    % SST
    data = data_format(filedir(j),'sst_d',row,lat,lon,3); %aicen
    w = worldmap('world');
        axesm eqaazim; 
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabelparallel', 0);
        setm(w, 'grid', 'on');
        setm(w, 'frame', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        colorbar
        title("SST",'interpreter','latex','FontSize', 12)
     nexttile
     % hi
     data = data_format(filedir(j),'hi_d',row,lat,lon,3); %aicen
    w = worldmap('world');
        axesm eqaazim; 
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabelparallel', 0);
        setm(w, 'grid', 'on');
        setm(w, 'frame', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        colorbar
        caxis([0,2]);
        title("Ice thickness",'interpreter','latex','FontSize', 12)
end

% afsdn


% for i = 1:24
%     idx = data(:,:,i)> 0;
%     sum(sum(idx))
% end
%data(lon_out,lat_out) = 10;
    
% for i = 1:4
%     [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
% end
if fstd_switch == 1
    data = data_format(filedir,'aicen_d',row,lat,lon,4);
    [lat_out,lon_out] = lat_lon_finder(coords(2,1),coords(2,2),lat,lon);
    fstd = data(lon_out,lat_out,:,:);
    size(sum(fstd))
end
    
    
    
    
%% Functions

