close all
clear all
addpath functions
% create video writer object
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'momentum';
grid = 'gx1'; 
variable = 'aice'; 
day = 3;
month = 7;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end
SIC = 0.15; 
dim = 2;
[lat,lon,row] = grid_read(grid);
filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
data_format(filename,"aice",row,lat,lon,dim);

obslat = 19.5362;
obslon = -155.5763;
obsH = 3397.00;

mllat = 19.475;
mllon = -155.608;
mlH = 4169;

figure
geoaxes("Basemap","satellite","ZoomLevel",12);
hold("on")
geoplot(obslat,obslon,"ow","MarkerSize",10,"MarkerFaceColor","magenta", ...
    "DisplayName","Mauna Loa Observatory");
geoplot(mllat,mllon,"ow","MarkerSize",10,"MarkerFaceColor","blue", ...
    "DisplayName","Mauna Loa Volcanao");
legend

figpos = [1000 500 800 400];
uif = uifigure("Position",figpos);
ug = uigridlayout(uif,[1,2]);
p1 = uipanel(ug);
p2 = uipanel(ug);
gx = geoaxes(p1,"Basemap","satellite"); 
gg = geoglobe(p2); 
gx.InnerPosition = gx.OuterPosition;
gg.Position = [0 0 1 1];


% View in 2-D
heightAboveTerrain = 200;
gx.MapCenter = [obslat, obslon];
zoomLevel = heightToZoomLevel(heightAboveTerrain, obslat);
gx.ZoomLevel = zoomLevel;

% View in 3-D
N = egm96geoid(obslat, obslon);
obsh = obsH + N;
ellipsoidalHeight = obsh + heightAboveTerrain;
campos(gg,obslat,obslon,ellipsoidalHeight)
drawnow


function zoomLevel = heightToZoomLevel(height, lat)
    earthCircumference = 2 * pi * 6378137;
    zoomLevel = log2((earthCircumference *cosd(lat)) / height) + 1;
    zoomLevel = max(0, zoomLevel);
    zoomLevel = min(19, zoomLevel);
end