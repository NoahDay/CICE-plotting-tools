clear all
close all
% Plotting the FSD hsitograms along a specified transect of the CICE output
% Comparing CICE results with and without waves
clear all
close all
filename = 'cases/1yearnowaves/history/iceh.2005-08-31.nc';
variable = "afsd";
latit = 1;
% Read the header
ncdisp(filename)
% grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');

dim = 2;
lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);


lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

% sea ice area
data_read = ncread(filename, variable); %dafsd_weld %dafsd_wave, dafsd_latm, dafsd_latg, dafsd_newi
fsd_cat = ncread(filename, 'NFSD');
fsd_range = [0; fsd_cat];
nfsd = length(fsd_cat);
%% Figure 1: for the first half of FSD categories.
figure(1)
%t1 = tiledlayout%(6,1);
for i = 1%:2%nfsd/2 % number of FSD categories
    data_2 = data_read(:,:,i);
    data = rearrange_matrix(data_2,37,dim);
    data = [data; data(end,:)];
    color_map = seaicecolormap();
    if max(max(data))>eps
        % Mapping
        %nexttile%(i)
        latitude = [-90,-30];
        longitude = [-180,180];
        w = worldmap('world');
        axesm eqdcylin; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
        setm(w, 'Origin', [0 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -85);
        setm(w, 'grid', 'on');
        %setm(w, 'frame', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        colorbar
        single_title = sprintf('[%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
        double_title = sprintf('FSD categories: \n [%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
        if i == 1
            title(double_title, 'Units', 'normalized', 'Position', [-0.15, 0.5, 0])
        else
            title(single_title, 'Units', 'normalized', 'Position', [-0.15, 0.5, 0]); % Set Title with correct Position
        end
        set(gca,'FontSize',10)
        %caxis([-0.01,0.01]);
    end
end

% Add layout title
if variable == "dafsd_wave"
 plot_var = " waves";
elseif variable == "dafsd_weld"
 plot_var = "welding";
end
plot_title = strcat('Change in the FSD due to ',plot_var, ' (FSD range [1,6])');
title(t1,plot_title)

%% Figure 2: for the second half of FSD categories.
figure(2)
t2 = tiledlayout(6,1);
for i = nfsd/2+1:nfsd % number of FSD categories
    data_2 = data_read(:,:,i);
    data = rearrange_matrix(data_2,37,dim);
    data = [data; data(end,:)];
    color_map = seaicecolormap();
    if max(max(data))>eps
        % Mapping
        nexttile%(i)
        latitude = [-90,-30];
        longitude = [-180,180];
        w = worldmap('world');
        axesm eqdcylin; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
        setm(w, 'Origin', [0 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -85);
        setm(w, 'grid', 'on');
        %setm(w, 'frame', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        colorbar
        single_title = sprintf('[%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
        double_title = sprintf('FSD categories: \n [%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
        if i == nfsd/2+1
            title(double_title, 'Units', 'normalized', 'Position', [-0.15, 0.5, 0])
        else
            title(single_title, 'Units', 'normalized', 'Position', [-0.15, 0.5, 0]); % Set Title with correct Position
        end
        set(gca,'FontSize',10)
    end
end

 % Add layout title
 if variable == "dafsd_wave"
     plot_var = "waves";
 elseif variable == "dafsd_weld"
     plot_var = "welding";
 end
 plot_title = strcat('Change in the FSD due to',plot_var,' (FSD range [7,12])');
 title(t2,plot_title)
