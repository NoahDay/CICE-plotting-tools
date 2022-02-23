function y = map_creator(filename, plot_title_vec, i, variable, grid, sector, user, SIC)
%MAP_CREATOR creates map images of Antarctica given netcdf data files
%   filename: the directory to the .nc file
%   plot_title: string containing the title for the plot
%   i: the date of the data file
%   variable: the variable in the dataset we want to plot
%   grid: specify the grid. eg. 'gx1', 'gx3'

% Load ice shelf data
addpath packages/bedmap2_toolbox_v4
%addpath ~/Users/noahday/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Bedmap2 Toolbox for Matlab


%shelf = bedmap2_data('icemask');
%shelf = 0.0*(shelf == 1); 
%shelf(shelf==0) = NaN; % 1 iceshelf
%[latshelf,lonshelf] = bedmap2_data('latlon');

% grid
    dim = 3;
[lat,lon,row] = grid_read(grid);

%filename = strcat('cases/',filename);
if variable == "vvel"
    data_x = ncread(filename, variable);
    data_y = ncread(filename, "uvel");
    data = sqrt(data_x.^2 + data_y.^2);
elseif variable == "vatm"
    data_x = ncread(filename, "uatm");
    data_y = ncread(filename, "vatm");
    data = sqrt(data_x.^2 + data_y.^2);
else
    if dim == 3
        data3d = ncread(filename, variable);
        data(:,:) = data3d(:,:,1);
    else
        data = ncread(filename, variable);
    end
end
[~, ~, n] = size(data);

    %% Mapping
    %color_map = seaicecolormap();
if sector == "world"
    x_origin = -90;
end
if sector == "world"
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [x_origin 0 0]);
        setm(w, 'maplatlimit', [-90,-30]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 30);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        %setm(w, 'frame', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        %pcolorm(latshelf,lonshelf,shelf)  
        colormap("cool")
elseif sector == "SA"
    [w a] = map_plot(data,variable,sector,grid);
    if SIC > eps
         aice = ncread(filename, "aice");
        [lat_ice_edge, lon_ice_edge] = find_ice_edge(aice,SIC,sector,lat,lon);
        plotm(lat_ice_edge,lon_ice_edge,'m-','LineWidth',2)
    end
end
 if variable == "fsdrad"
         plot_variable = "FSD radius ";
         unit = "metres";
     elseif  variable == "fsdrad_d"
         plot_variable = "FSD radius ";
         unit = "metres";
     elseif variable == "wave_sig_ht"
         plot_variable = "Significant wave height ";
         unit = "metres";
     elseif variable == "wave_sig_ht_d"
         plot_variable = "Significant wave height ";
         unit = "metres";
     elseif variable == "peak_period"
         plot_variable = "Peak period ";
         unit = "s";
     elseif variable == "peak_period_d"
         plot_variable = "Peak period ";
         unit = "s";
     elseif variable == "aice"
         plot_variable = "Concentration of ice ";
         unit = " ";
     elseif variable == "aice_d"
         plot_variable = "Concentration of ice ";
         unit = " ";
     elseif variable == "mean_wave_dir"
         plot_variable = "Mean wave direction (rads) ";
         unit = "radians";
     elseif variable == "mean_wave_dir_d"
         plot_variable = "Mean wave direction (rads) ";
         unit = "radians";
    elseif variable == "vvel"
        plot_variable = "Ice drift magnitude ";
        unit = "m/s";
    elseif variable == "vatm"
        plot_variable = "Wind magnitude ";
        unit = "m/s";
     else
         plot_variable = variable;
 end
    set(gcf, 'Position',  [100, 100, 1000, 800])
    set(gcf,'Visible', 'off')
    fontSize = 20; 
    plot_title = strcat(plot_variable, plot_title_vec);
    title(plot_title, 'FontSize', fontSize);
    %caxis([0 400]) %fsdrad [80 250]
    a=colorbar;
%    label_c = ylabel(a,unit,'FontSize',16,'Rotation',270);
    label_c.Position(1) = 4;
    label_h.Position(2) = 1; % change vertical position of ylabel
    limit = colorlims(variable);
    
    %labels = {'Forest','Water','Agriculture','Green Areas','Built-up'};
%lcolorbar(labels,'fontweight','normal', 'fontsize',16);
    caxis(limit);

    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %fname = '/Volumes/SSD/MATLAB/PhD Project/CICE Plotting/frames'; 
    %'/Users/noahday/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames';
    saveas(gcf,fullfile(filedir, figname));
end

