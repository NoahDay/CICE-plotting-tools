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
