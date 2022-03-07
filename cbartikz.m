%% Set custom colorbar ticks

%% Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar
clc
% Parameters
sector = "SA";
grid = 'gx1';
case_name = 'ocntest';
filedir = 'cases/ocntest/history/iceh.';%'/Volumes/NoahDay5TB/cases/ocnforcing/history/iceh.';
%'cases/momentum/history/iceh.'; %'/Volumes/NoahDay5TB/cases/momentum/history/iceh.2009-09-30.nc';
[lat,lon,row] = grid_read(grid);
user = "noahday";
% Make sector mask
%[len,wid] = size(lat);
fig_count = 0;


coords = sector_coords(sector);
% Take the location to be the centre of the sector
sector = "SA";%struct('coords',sector_coords(sector),'centre_lon',(coords(1,1)+coords(2,1))/2,'north_lat',[],'south_lat',[]);
clear coords

initial_date.day = 2;
initial_date.month = 1;
initial_date.year = 2009; % 2005 is a spin-up year
initial_date.char = sprintf('%d-0%d-0%d', initial_date.year, initial_date.month, initial_date.day);

filename = strcat(filedir,initial_date.char,'.nc');

pram.max_wind = 10; % m/s
pram.min_wind = 5; % m/s
pram.min_SIC = 0.15; % m/s
pram.total_datapoints = 1;
pram.quiver_color = 'red';
pram.label_color = 0.7*[0.4660 0.6740 0.1880];
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
pram.colormap = 'ice';


%fontname = 'Helvetica ';
%set(0,'defaultaxesfontname',fontname);
%set(groot,'defaulttextinterpreter','latex');  
%set(groot, 'defaultAxesTickLabelInterpreter','latex');  
%set(groot, 'defaultLegendInterpreter','latex'); 

aice_data = data_format_sector(filename,"aice",sector);
idx = aice_data < eps;
aice_data(idx) = NaN;
lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);

SIC = 0.15;
[lat_ice_edge, lon_ice_edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);

scale = 0;
 %colormap(summer(10))
%     cbh = colorbar('v');
%     get(get(cbh, 'Children'))
%     AxesH = axes('CLim', [0, 1000]);
%     %caxis([0 3000])
%     cbh = colorbar('peer', AxesH, 'v', ...
%                    'XTickLabel',{'0','100','1000'}, ...
%                    'XTick', [0,100,1000])
    %set(gca,'ColorScale','log')
    %colormap(turbo(20))
%    a = colorbar;
%            a.TickLabelInterpreter = 'latex';
%            a.Label.String = colorlabel(variable);
           % a.Ruler.Scale = 'log';
           

NFSD = ncread(filename,"NFSD");
c_tick = round([0;NFSD]);
% FSD radius
fsd_data = data_format_sector(filename,"fsdrad",sector);
idx = fsd_data < eps;
fsd_data(idx) = NaN;
figure(1)
    coords = sector_coords(sector);
        min_lat = min(coords(:,1));
        max_lat = max(coords(:,1));
        min_lon = min(coords(:,2));
        max_lon = max(coords(:,2));
        w = worldmap('world');
            axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            setm(w, 'Origin', [0 0 0]);
            setm(w, 'maplatlimit', [min_lat,max_lat]);
            setm(w, 'maplonlimit', [min_lon,max_lon]);
            setm(w, 'meridianlabel', 'on')
            setm(w, 'parallellabel', 'on')
            setm(w, 'mlabellocation', 10);
            setm(w, 'plabellocation', 10);
            setm(w, 'mlabelparallel', 0);
            setm(w, 'mlabelParallel', 'south');
            setm(w, 'grid', 'on');
            setm(w, 'frame', 'on');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,fsd_data)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.7 0.5]);
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1500 1000 500 600])
    cb = contourcbar('peer', w, 'v', ...
                       'XTickLabel',num2cell(c_tick'), ...
                       'XTick', c_tick, ...
                       'Location','eastoutside');
    cb.Ruler.Scale = 'log';
    colormap(turbo(10))
   get(get(cb, 'Children'))
   % AxesH = axes('CLim', [0, 1000]);
    %caxis([0 3000])
   % cbh = colorbar('peer', AxesH, 'v', ...
   %                'XTickLabel',{'0','100','1000'}, ...
   %                'XTick', [0,100,1000])



   load topo60c

   worldmap('north america')
contourfm(topo60c,topo60cR,-7000:1000:3000)
caxis([-8000 4000])
contourcbar
