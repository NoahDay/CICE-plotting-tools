clear all
close all
addpath functions
cd .. % Move up one directory
filename = 'cases/init/history/iceh_inst.2005-01-01-03600.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
latit = 1;
% Read the header
ncdisp(filename)
% 
% grid
lat = ncread('grid/global_gx1.bathy.nc','TLAT');
lon = ncread('grid/global_gx1.bathy.nc','TLON');

dim = 2;
lat = rearrange_matrix(lat,37,dim);
lon = rearrange_matrix(lon,37,dim);


lon = [zeros(1,384);lon];
lat = [lat(1,:); lat];

% sea ice area
% data = ncread(filename, 'afsd'); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
% data = data;
% %data = data(:,:,1);
% data = rearrange_matrix(data,37,dim);
% data = [data; data(end,:)];
% idx = data == 0;
% data(idx) = NaN;
%data(:,latit) = 10;
grid = 'gx1';
[lat,lon,row] = grid_read(grid);
data1 = data_format(filename,'afsd',row,lat,lon,3); %aicen

nfsd_cat = ncread(filename, 'NFSD');

for i = 1:320
    for j = 1:384
        temp = 0;
        for k = 1:24
            temp = temp + data1(i,j,k).*nfsd_cat(k);
        end
        if ~isnan(temp)
            if sum(temp) > eps
                data(i,j) = sum(temp);
            else
                 data(i,j) = NaN;
            end
        else
            data(i,j) = NaN;
        end
    end
end

%caxis([0,2])
color_map = seaicecolormap();
%% Mapping



latitude = [-90,-30];
longitude = [-180,180];


%idx = si_area_01 == 0;
%si_area_01 = si_area_01 + 0*idx;

%si_area_01 = rearrange_matrix(si_area_01,11);


%deglon = [deglon; deglon(100,:)+3.6];
%deglat = [deglat; deglat(100,:)];
%si_area_01 = [si_area_01; si_area_01(100,:)];


w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
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
    pcolorm(lat,lon,1000*data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

%antarctica = shaperead('landareas', 'UseGeoCoords', true,...
  %'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});


%load coastlines
%whos
%[latcells, loncells] = polysplit(coastlat, coastlon);
%numel(latcells)


%lon = linspace(0,360,100);
%lon = repmat(lon,116,1)';
%plotm(coastlat, coastlon)
%bedmap2('gl','xy')
%bedmap2('shelves')
%pcolorm(lat,lon,data)
%quiverm(lat,lon,data)
lat0 = [-90 -39.7]; 
lon0 = [-5.4 2.9];
deltalat = [5 5]; 
deltalon = [3 3];
%quiverm(lat0,lon0,deltalat,deltalon,'r')
colorbar
caxis([0,2000]);
%title("Representative FSD per cell (m) on 31-08-2005",'interpreter','latex','FontSize', 18)
%caxis([-1 1])
%patchm(antarctica.Lat, antarctica.Lon, [1 1 1])
%saveas(gcf,'test.png')