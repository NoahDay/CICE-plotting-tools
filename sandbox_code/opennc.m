clear all
close all
addpath functions
%cd .. % Move up one directory
day = 9;
month = 1;
year = 2005;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end
case_name = '31freq';
ticker = 1;
SIC = 0.15; 
filename = "/Users/noahday/GitHub/CICE-plotting-tools/cases/pancake_tracer/history/iceh.2005-06.nc";

%"/Users/noahday/GitHub/cice-dirs/input/CICE_data/ic/access-om2_1deg/iced.2019-01-01-00000.nc"%"/Users/noahday/GitHub/cice-dirs/runs/test/history/iceh.2005-01-24.nc";
%"/Users/noahday/github/cice-plotting-tools/cases/31freq/history/iceh.2008-07-30.nc";%/Volumes/NoahDay5TB/cases/forcingnowaves/history/iceh.2005-07.nc";%"prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_196201010130-196212312230.nc";
%strcat('cases/',case_name,"/history/iceh.",date,".nc");
% '/Users/noahday/Gadi/2010/JRA55_03hr_forcing_2010.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
%filename = 'DATA/CESM/MONTHLY/ocean_forcing_clim_2D_gx1.20210330.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';

% SWH
close all
conFigure(11,1.5)
lat = ncread(filename,"TLAT");
lon = ncread(filename,"TLON");
NFSD = ncread(filename,"NFSD");
data(:,:) = ncread(filename,"pancake_ice");
data4d(:,:,:,:) = ncread(filename,"afsdn");
data3d(:,:,:) = ncread(filename,"afsd");
data2(:,:) = data4d(:,:,1,1).*NFSD(1) ;
data3(:,:) = data3d(:,:,1);
dafsd_newi(:,:,:) = ncread(filename,"dafsd_newi");
%lat = ncread(filename,"TLAT");
%lon = ncread(filename,"TLON");

%data = data_format(filename,"wave_sig_ht");
%idx = data < eps;
%data(idx) = NaN;
figure;
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
    pcolorm(lat,lon,data(:,:))
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    caxis([0,1])
    cmocean('thermal',30)
%exportgraphics(f,'swhcice.pdf','ContentType','vector')
figure;
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
    pcolorm(lat,lon,data-data2)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    caxis([-1,1])
    cmocean('balance',11)

figure;
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
    pcolorm(lat,lon,data-data3)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    caxis([-1,1])
    cmocean('balance',11)

figure;
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
    pcolorm(lat,lon,dafsd_newi(:,:,1))
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    caxis([-0.01,0.01])
    cmocean('balance',11)

%% Ice age
close all
conFigure(11,1.5)
filename1= "cases/31freq/history/iceh.2008-07-31.nc";
filename2 = "cases/forcingoff/history/iceh.2008-07-31.nc";
f = figure
data1 = data_format(filename1,"iage");
data2 = data_format(filename2,"iage");
%idx = data1 < 12/12;
%data1(idx) = NaN;
[p,a] = map_plot(data1-data2,"iage","SH","gx1",[-0.5,0.5])
cmocean('balance')
a.Label.String = "Ice age, years";
a.Label.Interpreter = "latex";
%a.Limits = [-0.5,0.5];

exportgraphics(f,'iagecomp.pdf','ContentType','vector')

%%
f = figure;
idx = data < 12/12;
data(idx) = NaN;
[p,a] = map_plot(data,"iage","SH")
a.Label.String = "Ice age, years";
a.Label.Interpreter = "latex";
%title("Ice older than 6 months")
%exportgraphics(f,'iage12monthnoforc.pdf','ContentType','vector')
%% SWH
filename = "/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/iceh.2006-07.nc";
close all
conFigure(11,1.5)
f = figure
data = data_format(filename,"wave_sig_ht");
idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"wave_sig_ht","SH")
s = scaleruler
setm(handlem('scaleruler'), ...
    'XLoc',0.3, ... '
    'YLoc',-0.52, ...
    'TickDir','down', ...
    'MajorTick',0:1000:1000, ...
    'MinorTick',0:500:500, ...
    'MajorTickLength',km2nm(150),...
    'MinorTickLength',km2nm(150))
a.Label.String = "Signficant wave height (m)";
a.Label.Interpreter = "latex";
%exportgraphics(f,'swhcice.pdf','ContentType','vector')

idx = isnan(data);
mean(data(~idx))
%%
f = figure
data = data_format(filename,"peak_period");
idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"peak_period","SH")
s = scaleruler
setm(handlem('scaleruler'), ...
    'XLoc',0.3, ... '
    'YLoc',-0.52, ...
    'TickDir','down', ...
    'MajorTick',0:1000:1000, ...
    'MinorTick',0:500:500, ...
    'MajorTickLength',km2nm(150),...
    'MinorTickLength',km2nm(150))
a.Label.String = "Peak period (s)";
a.Label.Interpreter = "latex";
%exportgraphics(f,'ppdcice.pdf','ContentType','vector')

idx = isnan(data);
mean(data(~idx))
%% Concentration of small floes
close all
clear data3d
clear data3d1 data3d2
f = figure;

NCAT = ncread(filename,"NCAT");
NFSD = ncread(filename,"NFSD");
Nf = numel(NFSD);

floe_binwidth = [5.2438,8.9763,14.7711,23.3545,35.4569,51.6493,72.1173,96.4015,123.1658,150.0742,173.8638,190.6718,397.7316,479.1093,649.9598,881.7363];
thick_binwidth = NCAT - [0;NCAT(1:end-1)];


data3d1(:,:,:) = data_format(filename1,"afsd");
data3d2(:,:,:) = data_format(filename2,"afsd");
data(:,:) = data3d1(:,:,1).*floe_binwidth(1)-data3d2(:,:,1).*floe_binwidth(1);
idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"aice","SH","gx1",[-1,1]);
% s = scaleruler
% setm(handlem('scaleruler'), ...
%     'XLoc',0.3, ... '
%     'YLoc',-0.52, ...
%     'TickDir','down', ...
%     'MajorTick',0:1000:1000, ...
%     'MinorTick',0:500:500, ...
%     'MajorTickLength',km2nm(150),...
%     'MinorTickLength',km2nm(150))
a.Label.String = "Ice concentration";
a.Label.Interpreter = "latex";
cmocean('balance',10)
exportgraphics(f,'pancakeconccomp.pdf','ContentType','vector')

clear data3d1 data3d2
f = figure;
data3d1 = data_format(filename1,"aicen");
data3d2 = data_format(filename2,"aicen");

data(:,:) = data3d1(:,:,1)-data3d2(:,:,1);

idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"aice","SH","gx1",[-1,1])
% s = scaleruler
% setm(handlem('scaleruler'), ...
%     'XLoc',0.3, ... '
%     'YLoc',-0.52, ...
%     'TickDir','down', ...
%     'MajorTick',0:1000:1000, ...
%     'MinorTick',0:500:500, ...
%     'MajorTickLength',km2nm(150),...
%     'MinorTickLength',km2nm(150))
a.Label.String = "Ice concentration";
a.Label.Interpreter = "latex";
cmocean('balance',10)
exportgraphics(f,'thinconccomp.pdf','ContentType','vector')

clear data4d1 data4d2
f = figure;
data4d1 = data_format(filename1,"afsdn");
data4d2 = data_format(filename2,"afsdn");
% for i = 1:321
%     for j = 1:384
%         temp = data3d(i,j,1);
%         data(i,j) = 
data(:,:) = data4d1(:,:,1,1).*floe_binwidth(1)-data4d2(:,:,1,1).*floe_binwidth(1);
idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"aice","SH","gx1",[-1,1])
% s = scaleruler
% setm(handlem('scaleruler'), ...
%     'XLoc',0.3, ... '
%     'YLoc',-0.52, ...
%     'TickDir','down', ...
%     'MajorTick',0:1000:1000, ...
%     'MinorTick',0:500:500, ...
%     'MajorTickLength',km2nm(150),...
%     'MinorTickLength',km2nm(150))
a.Label.String = "Ice concentration";
a.Label.Interpreter = "latex";
cmocean('balance',10)
exportgraphics(f,'smallthinconccomp.pdf','ContentType','vector')

%%
close all
clear data p
filename_x = "cases/forcingoff/history/iceh.2005-07-31.nc";
filename_x = "/Volumes/NoahDay5TB/cases/noforcing/history/iceh.2005-01.nc";
%filename_x = "cases/forcingoff/history/iceh_ic.2005-01-01-03600.nc";
data(:,:) = data_format(filename_x,"fsdrad");
idx = data < eps;
data(idx) = NaN;
[p,a] = map_plot(data,"fsdrad","SH","gx1",[0,1000])
% s = scaleruler
% setm(handlem('scaleruler'), ...
%     'XLoc',0.3, ... '
%     'YLoc',-0.52, ...
%     'TickDir','down', ...
%     'MajorTick',0:1000:1000, ...
%     'MinorTick',0:500:500, ...
%     'MajorTickLength',km2nm(150),...
%     'MinorTickLength',km2nm(150))
%p.ZScale = 'log';
a.Label.String = "Floe size radius [m]";
a.Label.Interpreter = "latex";
cmocean('ice',10)
%exportgraphics(f,'fsdradnoforc.pdf','ContentType','vector')

%%
%cd /Users/a1724548/Github/CICE-plotting-tools
f3 = figure
ppddata = data_format(filename,"peak_period");
[p,a] = map_plot(ppddata,"peak_period","SH")
a.Label.String = "Peak period, $s$";
a.Label.Interpreter = "latex";
 %exportgraphics(f3,'peakperiod.pdf','ContentType','vector')
colormap parula

idx = ppddata > eps;
f = figure;
hist(ppddata(idx))
xlabel('Peak period, $s$','Interpreter','latex')
ylabel('Count','Interpreter','latex')
%exportgraphics(f,'hist.pdf','ContentType','vector')
% SWH
f3 = figure;
swhdata = data_format(filename,"wave_sig_ht");
[p,a] = map_plot(swhdata,"wave_sig_ht","SH");
a.Label.String = "SWH, $m$";
a.Label.Interpreter = "latex";
%exportgraphics(f3,'swh.pdf','ContentType','vector')
colormap parula

idx = swhdata > eps;
f = figure;
hist(swhdata(idx))
xlabel('Significant wave height, $m$','Interpreter','latex')
ylabel('Count','Interpreter','latex')
%exportgraphics(f,'histswh.pdf','ContentType','vector')


%% 
close all
conFigure(11,1.2)
nbins = 50;
f = figure;
historydir = '/Users/noahday/GitHub/CICE-plotting-tools/cases/31freq/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

a = dir([historydir '/*.nc']);
n_files = numel(a);
swhdata_vec = [];
ppddata_vec = [];
for i = 1:n_files
    clear swhtemp ppdtemp
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
   swhdata_mat(:,:,i) = data_format(filenames(i,:),"wave_sig_ht");
   ppddata_mat(:,:,i) = data_format(filenames(i,:),"peak_period");
   swhtemp(:,:) = swhdata_mat(:,:,i);
   swhdata_vec = [swhdata_vec; swhtemp(idx)];
   ppdtemp(:,:) = ppddata_mat(:,:,i);
   ppddata_vec = [ppddata_vec; ppdtemp(idx)];
end




hist3([swhdata_vec, ppddata_vec],'Nbins',[nbins,nbins],'CdataMode','auto')

xlabel('Significant wave height, $H_s$ (m)','Interpreter','latex')
    xlim([0,max(swhdata(idx))])
    ylim([0,max(ppddata(idx))])
    ylabel('Peak period, $T_p$ (s)','Interpreter','latex')
    zlabel('Counts','Interpreter','latex')
    s.EdgeColor = 'none';
    a = colorbar;
    a.Label.String = 'Counts';
    a.Label.Interpreter = 'latex';
    a.TickLabelInterpreter = 'latex';
    C=jet;
    C(1,:) = [1,1,1];
    colormap(C)
    caxis([eps 1200])
    view(0,90)
    
exportgraphics(f,'countswhppd.pdf','ContentType','vector')
%%
surf(swhdata(idx), ppddata(idx))

%% change in new ice
close all
clear f
%conFigure(50)
filename = "cases/forcingoff/history/iceh.2008-07-31.nc";
f = figure;
dafsd3d = data_format(filename,"dafsd_");
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
for i = 1:16
    nexttile
    dafsd = dafsd3d(:,:,i)*NFSD(i);%.*floe_binwidth(i);
    idx = dafsd < eps;
    dafsd(idx) = NaN;
    if isnan(max(max(dafsd)))
        maxi = 1;
    else
        maxi = max(max(dafsd));
    end
    [p,a] = map_plot(dafsd,"dafsd_wave","SH","gx1",[0,maxi]);
    if isnan(max(max(dafsd)))
        title('No Change','FontSize',12)
    else
        maxi = max(max(dafsd));
    end
    %a.Label.String = "Change in FSD";
    a.Label.Interpreter = "latex";
    colormap parula
end
t.YLabel.String = "Change in areal concentration (dimensionless)";
t.YLabel.Interpreter = "latex";
t.Ylabel.FontSize = a.FontSize;
% exportgraphics(f,'dafsd.pdf','ContentType','vector')
%%
idx = fsddata > eps;
f = figure;
hist(fsddata(idx))
xlabel('Significant wave height, $m$','Interpreter','latex')
ylabel('Count','Interpreter','latex')
exportgraphics(f,'histswh.pdf','ContentType','vector')
%%
%lat = ncread('grid/gx1/global_gx1.bathy.nc','TLAT');
%lon = ncread('grid/gx1/global_gx1.bathy.nc','TLON');
%[lat,lon,row] = grid_read(grid);
%data = ncread(filename, 'U');
%data = data_format(filename,'prra',row,lat,lon,3);
%data = data(:,:,1);
[lat,lon,row] = grid_read(grid);
lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);
dim = 2;

data_u = data_format(filename,'strocnx',row,lat,lon,dim);
data_u2 = data_u(:,:,1);
data_v = data_format(filename,'strocny',row,lat,lon,dim);
data_v2 = data_v(:,:,1);
u_vec = reshape(data_u2,1,[]);
v_vec = reshape(data_v2,1,[]);

data = sqrt(data_u2.^2 + data_v2.^2);
data =  data_format(filename,'sst',row,lat,lon,dim);
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
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %quiverm(lat_vec,lon_vec,u_vec,v_vec,'k')
    
  %%
close all
filename = '/Users/noahday/Maths1/iceh.2005-01-01.nc';

lat = ncread(filename,"ULAT");
lon = ncread(filename,"ULAT"); 
[lat,lon,row] = grid_read("gx1");
swh = data_format_sector(filename,"wave_sig_ht","SH");
aice = data_format(filename,"aice");
idx = aice > 0.01;
swh(~idx) = 0;
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
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,swh)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('balance',31)
    caxis([-3,3])


    %%
    
