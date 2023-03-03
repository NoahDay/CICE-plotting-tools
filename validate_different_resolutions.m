clear all
close all


filename_10 = '/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh.2017-01-01.nc';

lat_10 = ncread(filename_10,'TLAT');
lon_10 = ncread(filename_10,'TLON');

tair_10 = ncread(filename_10,'Tair');


filename_025 = '/Users/noahday/GitHub/cice-dirs/runs/access-om2_025/history/iceh_inst.2003-01-01-03600.nc';

lat_025 = ncread(filename_025,'TLAT');
lon_025 = ncread(filename_025,'TLON');

tair_025 = ncread(filename_025,'Tair');


[len wid] = size(lat_10);
for i = 1:len
    for j = 1:wid
        lat_10_025(i,1) = mean(lat_025(4*(i-1)+1:4*(i-1)+4,1));
    end
end


conFigure(11,1)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
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
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat_10,lon_10,tair_10) 
    a = colorbar;
    a.Label.String = "Air temp [C]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([-30,30])


conFigure(11,1)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
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
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat_025,lon_025,tair_025) 
    a = colorbar;
    a.Label.String = "Air temp [C]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([-30,30])


%% Check the the ocean velocities were read in correctly on Maths1
% We simply check whether the forcing files change each day

filename1 = '/Users/noahday/GitHub/cice-dirs/runs/access-om2_025/history/iceh_inst.2003-01-01-03600.nc';
filename2 = '/Users/noahday/GitHub/cice-dirs/runs/access-om2_025/history/iceh_inst.2003-01-01-07200.nc';


[lat,lon] = grid_read('om2');
lat = ncread(filename1,'ULAT');
lon = ncread(filename1,'ULON');
uocn1 = ncread(filename1,'uocn');
uocn2 = ncread(filename2,'uocn');


conFigure(11,1)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
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
    setm(w, 'grid', 'on');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,uocn2-uocn1) 
    a = colorbar;
    a.Label.String = "Difference in u-ocean velocity [m/s]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([-0.001,0.001])
    cmocean('balance')
exportgraphics(f,'difference_u-ocean.png','ContentType','image')
