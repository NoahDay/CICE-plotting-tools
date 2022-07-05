%% Case info
clear all
clc
close all
addpath plotting
addpath processing

historydir = '/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end
lat = data_format(filenames(1,:),"TLAT");

sector = "SH";
for i = 1:1 %n_files
    afsdn = data_format_sector(filenames(i,:),"afsdn",sector);
    afsd = data_format_sector(filenames(i,:),"afsd",sector);
    dafsd_latg = data_format_sector(filenames(i,:),"dafsd_latg",sector);
    dafsd_latm = data_format_sector(filenames(i,:),"dafsd_latm",sector);
    dafsd_newi = data_format_sector(filenames(i,:),"dafsd_newi",sector);
    dafsd_weld = data_format_sector(filenames(i,:),"dafsd_weld",sector);
    dafsd_wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
end

pancake_ice(:,:) = afsd(:,:,1) - dafsd_wave(:,:,1);

figure;
w = worldmap('world');
    axesm miller; 
    setm(w, 'Origin', [90 0 0]);
    setm(w, 'maplatlimit', [60,90]);
    setm(w, 'maplonlimit', [0,360]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'on');
    setm(w, 'labelrotation', 'on')
    plotm(grid.access.lat(1:10:end,1:10:end),grid.access.lon(1:10:end,1:10:end),10)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    title('ACCESS-OM2 grid')