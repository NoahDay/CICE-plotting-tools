%% FSD rad
close all
clear all
clc
addpath functions
filename_wim = "/Volumes/NoahDay5TB/cases/monthwim/history/iceh.2005-09.nc";
filename_none = "/Volumes/NoahDay5TB/cases/noforcing/history/iceh.2005-07.nc";
grid = "gx1";
[lat,lon,row] = grid_read(grid);

lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);
dim = 2;


SIC_threshold = 0.15;
sector = "SH";
aice_wim = data_format(filename_wim,'aice');
afsdn_none = data_format(filename_none,'afsdn');
aice_none = fsd_converter(filename_none,"afsdn","aice",afsdn_none,sector);

% Ice mask
icemask_wim = aice_wim > SIC_threshold;
icemask_none = aice_none > SIC_threshold;
%% Plotting
close all
clc
% Representative radius WITH waves
data = data_format(filename_wim,'fsdrad');
idx = data < eps;
data(idx) = NaN;
data(~icemask_wim) = NaN;
conFigure(11,1)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Mean floe size radius [m]';
    c.Label.Interpreter = 'latex';
    caxis([0,3000]);
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    %exportgraphics(f,'swhtest.pdf','ContentType','vector')
    
 % Representative radius WITHOUT waves   
data = data_format(filename_none,'fsdrad');
idx = data < eps;
data(idx) = NaN;
data(~icemask_none) = NaN;
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Mean floe size radius [m]';
    c.Label.Interpreter = 'latex';
    caxis([0,3000]);
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    
% Significant wave height WITH waves        
data = data_format(filename_wim,'wave_sig_ht');
idx = data < eps;
data(idx) = NaN;
data(~icemask_wim) = NaN;
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Mean floe size radius [m]';
    c.Label.Interpreter = 'latex';
    %c.Limits = [0,5];
    caxis([0,5])
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    
%% Change in FSD    
close all
% Waves      
dafsd_wave = data_format(filename_wim,'dafsd_wave');
data = fsd_converter(filename_wim,"dafsd_wave","fsdrad",dafsd_wave,sector);
data(~icemask_wim) = NaN;
data = data .*aice_wim;
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Change in FSD due to waves [m/day]';
    c.Label.Interpreter = 'latex';
    cmocean('balance')
    %c.Limits = [0,5];
    caxis([-10,10])
    
    
% New ice      
dafsd_wave = data_format(filename_wim,'dafsd_newi');
data = fsd_converter(filename_wim,"dafsd_newi","fsdrad",dafsd_wave,sector);
data(~icemask_wim) = NaN;
data = data .*aice_wim;
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Change in FSD due to new floes [m/day]';
    c.Label.Interpreter = 'latex';
    cmocean('balance')
    %c.Limits = [0,5];
    caxis([-10,10])
    
    
%%
%Tilt
datax = data_format(filename_wim,'strtltx');
datay = data_format(filename_wim,'strtlty');
data = sqrt(datax.^2 + datay.^2);
%data = fsd_converter(filename_wim,"dafsd_newi","fsdrad",dafsd_wave,sector);
data(~icemask_wim) = NaN;
data = data .*aice_wim;
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-55]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    %setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Sea surface tilt (N/m$^2$)';
    c.Label.Interpreter = 'latex';
    cmocean('balance')
    %c.Limits = [0,5];
    caxis([-0.0291,0.03])
    exportgraphics(f,'tilt.pdf','ContentType','vector')


%%
filename = '/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh.2017-01-01.nc';
[lat,lon] = grid_read('om2');

tmask = data_format(filename,'tmask');

lat = ncread(filename,'TLAT');
lon = ncread(filename,'TLON');
tmask = ncread(filename,'tmask');

set(gcf, 'color', 'none')

close all
f = figure;
latlim = [-90 0];
lonlim = [-90 90];
w = worldmap('world'); 
    %axesm eqaazim
    setm(w, 'Origin', [0 -180 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-270,90]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    %setm(w, 'mlabellocation', 60);
    %setm(w, 'plabellocation', 10);
    %setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    %setm(w, 'labelrotation', 'on')
    colormap([0.8, 0.8, 0.8; 1,1,1])
    pcolorm(lat,lon,tmask)
    %plotm(lat(:,1:10:end),lon(:,1:10:end),'b')
    %hold on
    for i = 1:300
        plotm(lat(1,:),lon(i,:),'b')
    end
    for i = 1:5:length(lat(:,1))
       % plotm(lat(:,i),lon(:,1),'b')
    end
    %plotm(lon,lat,'b')
    
%exportgraphics(f,'1deg_grid.pdf','ContentType','vector')

%%
close all
f = figure;
latlim = [50 59];
lonlim = [-10.5 2];
w = worldmap(latlim,lonlim); 
    %setm(w, 'Origin', [40 0 0]);
    %setm(w, 'maplatlimit', [-30,30]);
    %setm(w, 'maplonlimit', [50,70]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    %setm(w, 'mlabellocation', 60);
    %setm(w, 'plabellocation', 10);
    %setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    %setm(w, 'labelrotation', 'on')
    colormap([0.8, 0.8, 0.8; 1,1,1])
    pcolorm(lat,lon,tmask)
    %plotm(lat(:,1:10:end),lon(:,1:10:end),'b')
    %hold on
    %for i = 1:10:length(lat(1,:))
    %    plotm(lat(i,:),lon(i,:),'b')
    %end
    %plotm(lon,lat,'b')
    title('1$^\circ$ resolution','FontSize',20, 'Position', [0.5, 0.75, 0]);


%matlab2tikz('Britain_1deg.tex', 'standalone', true)
exportgraphics(f,'Britain_1deg.pdf','ContentType','vector')
%imwrite(bitmapData, 'Britain_1deg.png', 'png', 'transparency', backgroundColor)

%f = figure;


filename_025 = '/Users/noahday/GitHub/cice-dirs/runs/access-om2_025/history/iceh_inst.2003-01-01-03600.nc';

lat_025 = ncread(filename_025,'TLAT');
lon_025 = ncread(filename_025,'TLON');
tmask_025 = ncread(filename_025,'tmask');

f = figure;
w = worldmap(latlim,lonlim);
setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    %setm(w, 'labelrotation', 'on')
    colormap([0.8, 0.8, 0.8; 1,1,1])
    pcolorm(lat_025,lon_025,tmask_025)
    title('0.25$^\circ$ resolution','FontSize',20, 'Position', [0.5, 0.75, 0]);

%matlab2tikz('Britain_025deg.tex', 'standalone', true)
exportgraphics(f,'Britain_025deg.pdf','ContentType','vector')
%imwrite(bitmapData, 'Britain_025deg.png', 'png', 'transparency', backgroundColor)
% save to transparented image
%set(gcf, 'color', 'none');    
%set(gca, 'color', 'none');
%exportgraphics(gcf,'transparent.emf',...   % since R2020a
%    'ContentType','vector',...
%    'BackgroundColor','none')