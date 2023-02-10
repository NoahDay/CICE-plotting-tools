% Check that the CICE simulation is running correctly

clc
clear all
close all
addpath functions

filename = '/Users/noahday/GitHub/cice-dirs/runs/access-om2_025/history/iceh_inst.2003-01-02-00000.nc';


sector = 'SH';
grid = "om2_025";
%aice_data = data_format_sector(filename,var,sector);
%[lat,lon] = grid_read(grid);

lat = ncread(filename,'TLAT');
lon = ncread(filename,'TLON');

conFigure(11,1)
%%
var = 'aice';
aice_data = ncread(filename,var);
[f, a] = make_map(lat,lon,aice_data);
    cmocean('ice')
    a.Label.String = "Sea ice concentration [-]";
    caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear aice_data

%% Atmospheric forincg
set(gcf,'Visible', 'off');
var = 'Tair';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('thermal')
    a.Label.String = "Surface air temperature [$^\circ$ C]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data

var = 'uatm';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('balance','pivot',0)
    a.Label.String = "Zonal wind velocity [m/s]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data

var = 'vatm';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('balance','pivot',0)
    a.Label.String = "Meridional wind velocity [m/s]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data

var = 'snow';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('tempo')
    a.Label.String = "Snow fall rate [cm/day]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data

var = 'rain';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('tempo')
    a.Label.String = "Rainfall rate [cm/day]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data


%%
var = 'uvel';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('balance','pivot',0)
    a.Label.String = "Zonal ice velocity [m/s]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data

%% 
var = 'uvel';
data = ncread(filename,var);
[f, a] = make_map(lat,lon,data);
    cmocean('balance','pivot',0)
    a.Label.String = "Zonal ice velocity [m/s]";
    %caxis([0,1])
exportgraphics(f,strcat(var,'_map.png'),'ContentType','image',"Resolution",2000)
clear data



















% FUNCTIONS

function [f_out, a] = make_map(lat,lon,data)
    f_out = figure;
    w = worldmap('world');
        axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-50]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'off')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 60);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'mlinelimit', [-75 -55]);
        setm(w, 'plinelimit', [-75 -55]);
        setm(w, 'grid', 'off');
        setm(w, 'frame', 'off');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data) 
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        a = colorbar;
        a.TickLabelInterpreter = 'latex';
        a.Label.Interpreter = 'latex';
end