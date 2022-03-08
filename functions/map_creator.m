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
[lat,lon,row] = grid_read(grid);

dim = 2;
%filename = strcat('cases/',filename);
if variable == "vvel"
    dim = 2;
    data_x = data_format(filename,"uvel",row,lat,lon,dim);
    data_y = data_format(filename,"vvel",row,lat,lon,dim);
    data = sqrt(data_x.^2 + data_y.^2);
elseif variable == "vatm" 
    dim = 2;
    data_x = data_format(filename,"uatm",row,lat,lon,dim);
    data_y = data_format(filename,"vatm",row,lat,lon,dim);
    data = sqrt(data_x.^2 + data_y.^2);
elseif variable == "Tsfc" || variable == "Tair" || variable == "sst"
    dim = 2;
    data_x = data_format(filename,"uatm",row,lat,lon,dim);
    data_y = data_format(filename,"vatm",row,lat,lon,dim);
    data = data_format(filename,variable,row,lat,lon,dim);
elseif variable == "stresses"
    lat_vec = reshape(lat,1,[]);
    lon_vec = reshape(lon,1,[]);
    % Read in all the stress data and format into vectors
    strairx_data = data_format_sector(filename,"strairx",sector);
    strairy_data = data_format_sector(filename,"strairy",sector);
    
    strairx_vec = reshape(strairx_data,1,[]);
    strairy_vec = reshape(strairy_data,1,[]);
    
    strocnx_data = data_format_sector(filename,"strocnx",sector);
    strocny_data = data_format_sector(filename,"strocny",sector);
    
    strocnx_vec = reshape(strocnx_data,1,[]);
    strocny_vec = reshape(strocny_data,1,[]);
    
    strcorx_data = data_format_sector(filename,"strcorx",sector);
    strcory_data = data_format_sector(filename,"strcory",sector);
    
    strcorx_vec = reshape(strcorx_data,1,[]);
    strcory_vec = reshape(strcory_data,1,[]);
    
    
    strtltx_data = data_format_sector(filename,"strtltx",sector);
    strtlty_data = data_format_sector(filename,"strtlty",sector);
    
    strtltx_vec = reshape(strtltx_data,1,[]);
    strtlty_vec = reshape(strtlty_data,1,[]);

    strintx_data = data_format_sector(filename,"strintx",sector);
    strinty_data = data_format_sector(filename,"strinty",sector);

    strintx_vec = reshape(strintx_data,1,[]);
    strinty_vec = reshape(strinty_data,1,[]);

    data = data_format(filename,"aice",row,lat,lon,dim);
    idx = data < 0.01;
    data(idx) = NaN;
else
    if dim == 3
        data3d = ncread(filename, variable);
        data(:,:) = data3d(:,:,1);
    else
        data = data_format(filename,variable,row,lat,lon,dim);%ncread(filename, variable);
    end
end


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
    if variable == "fsdrad"
        data = data_format(filename,"fsdrad",row,lat,lon,2);
        idx = data < 0.01;
        data(idx) = NaN;% = data.*idx;
    end
    [w,a] = map_plot(data,variable,sector);
    
    if SIC > eps
        aice = data_format(filename,"aice",row,lat,lon,2);
        idx = aice < 0.01;
        aice(idx) = NaN;
        [lat_ice_edge, lon_ice_edge] = find_ice_edge(aice,SIC,sector,lat,lon);
        t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',SIC),'HorizontalAlignment','right');
    end
    if variable == "vatm" || variable == "Tair" || variable == "Tsfc" || variable == "sst"
        plotm(lat_ice_edge,lon_ice_edge,'-','color',[0.99 0.99 0.99],'LineWidth',2)
        lat_vec = reshape(lat,1,[]);
        lon_vec = reshape(lon,1,[]);
        x_vec = reshape(data_x,1,[]);
        y_vec = reshape(data_y,1,[]);
        color = 'thermal';
        cmocean(color,15)
        quiverm(lat_vec,lon_vec,x_vec,y_vec,'cyan') 
    elseif variable == "vvel"
        lat_vec = reshape(lat,1,[]);
        lon_vec = reshape(lon,1,[]);
        x_vec = reshape(data_x,1,[]);
        y_vec = reshape(data_y,1,[]);
        quiverm(lat_vec,lon_vec,x_vec,y_vec,'cyan') 
        plotm(lat_ice_edge,lon_ice_edge,'-','color',[0.99 0.99 0.99],'LineWidth',2)
        colormap(turbo(20))
    elseif variable == "wave_sig_ht"
        plotm(lat_ice_edge,lon_ice_edge,'-','color',[0.99 0.99 0.99],'LineWidth',2)
        colormap(winter(20))
    elseif variable == "fsdrad"
        plotm(lat_ice_edge,lon_ice_edge,'-','color',[0.99 0.99 0.99],'LineWidth',2)
        %colormap(jet(30))
    elseif variable == "aice"
        plotm(lat_ice_edge,lon_ice_edge,'-','color',0.7*[0.4660 0.6740 0.1880],'LineWidth',2) % Dark green
        cmocean('ice',10)
    elseif variable == "stresses"
        scale = 0;
        plotm(lat_ice_edge,lon_ice_edge,'-','color',0.7*[0.4660 0.6740 0.1880],'LineWidth',2)
        q1 = quiverm(lat_vec,lon_vec,strairx_vec,strairy_vec,'b',scale);
        q2 = quiverm(lat_vec,lon_vec,strocnx_vec,strocny_vec,'r',scale);
        q3 = quiverm(lat_vec,lon_vec,strcorx_vec,strcory_vec,'g',scale);
        q4 = quiverm(lat_vec,lon_vec,strtltx_vec,strtlty_vec,'m',scale);
        q5 = quiverm(lat_vec,lon_vec,strintx_vec,strinty_vec,'k',scale);
        legend([q1(1), q2(1), q3(1), q4(1), q5(1)], 'Air', 'Ocean','Coriolis','Sea surface slope','Internal','Box','on','NumColumns',2);
         %legend([q1(1), q2(1), q3(1), q4(1), q5(1)], 'Air', 'Ocean','Coriolis','Sea surface slope','Internal')
        %pos = get(h,'Position');
        posx = 0.35;
        posy = 0.75;
        cmocean('ice',10)  
        %set(h,'Position',[posx posy pos(3) pos(4)]);
    else
        plotm(lat_ice_edge,lon_ice_edge,'-','color',[0.99 0.99 0.99],'LineWidth',2)
        colormap(turbo(20))
    end
    
      %  edge_color = 0.7*[0.4660 0.6740 0.1880]; % Dark green
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
        plot_variable = "Wind forcing ";
        unit = "m/s";
    elseif variable == "sst"
        plot_variable = "Sea surface temperature ";
        unit = "C";
    elseif variable == "Tair"
        plot_variable = "Air temperature ";
        unit = "C";
 elseif variable == "stresses"
     plot_variable = "All ice stresses ";
     unit = "Concentration";
     else
         plot_variable = variable;
 end
    set(gcf, 'Position',  [100, 100, 400, 420])
    set(gcf,'Visible', 'off')
    fontSize = 16; 
    plot_title = strcat(plot_variable, plot_title_vec);
    title(plot_title, 'FontSize', fontSize);

    figname = sprintf('image%d.png', i); 
    filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
    %fname = '/Volumes/SSD/MATLAB/PhD Project/CICE Plotting/frames'; 
    %'/Users/noahday/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames';
    saveas(gcf,fullfile(filedir, figname));
end

