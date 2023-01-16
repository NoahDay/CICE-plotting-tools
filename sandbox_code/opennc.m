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
[lat,lon] = grid_read('om2');

conFigure(11,1.5)
filename1= "/Volumes/NoahDay5TB/WIMonAlessandroRun/history/iceh.2017-09-08.nc";
%filename2 = "/Volumes/NoahDay5TB/WIMoffSheetIceAlessandroRun/history/iceh.2019-12-31.nc";
f1 = figure;
data1 = data_format(filename1,"iage");
%data2 = data_format(filename2,"iage");
%idx = data1 < 12/12;
%data1(idx) = NaN;
%[p,a] = map_plot(data1-data2,"iage","SH","gx1",[-0.5,0.5])
% figure;
% w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [-90 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'on')
%     setm(w, 'parallellabel', 'off')
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     %setm(w, 'frame', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(lat,lon,data1-data2)
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     colorbar
%     caxis([-1,1])
%     cmocean('balance','pivot')
% a.Label.String = "Ice age, years";
% a.Label.Interpreter = "latex";
%a.Limits = [-0.5,0.5];

%exportgraphics(f,'iagecomp.pdf','ContentType','vector')



f2 = figure;
data1 = data_format(filename1,"aice");
%data2 = data_format(filename2,"aice");
%idx = data1 < 12/12;
%data1(idx) = NaN;
%[p,a] = map_plot(data1-data2,"iage","SH","gx1",[-0.5,0.5])
f = figure;
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
    pcolorm(lat,lon,data1)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    a = colorbar
    caxis([0.9,1])
    cmocean('ice')
a.Label.String = "SIC";
a.Label.Interpreter = "latex";

exportgraphics(f,'SIC90.pdf','ContentType','vector')
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
filename ='/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history_h/iceh_inst.2017-07-03-57600.nc'
close all
conFigure(11,1.5)
f = figure
data = data_format(filename,"Tair");
%idx = data < eps;
%data(idx) = NaN;
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
%caxis(0,3)
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
historydir = '/Users/noahday/Maths1/access-forcing-2010-new-ocean-settings/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

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
filename = "/Users/noahday/Maths1/access-forcing-2010-new-ocean-settings/iceh.2010-08.nc";
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
filename = '/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh_inst.2017-07-01-82800.nc';

lat = ncread(filename,"ULAT");
lon = ncread(filename,"ULON"); 
[latgx1,longx1,row] = grid_read("gx1");
swh = data_format_sector(filename,"wave_sig_ht","vichi"); %data_format_sector(filename,"wave_sig_ht","SH");
aice = data_format(filename,"aice");
idx = aice > 0.01;
%swh(~idx) = 0;
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
    caxis([-8,8])


    %%
close all
conFigure(30,3)
figure
subplot(1,3,1)
filename = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history_h/iceh_inst.2017-07-03-82800.nc'
[lat,lon,row] = grid_read("om2");
swh2 = data_format(filename,"frzmlt"); %data_format_sector(filename,"wave_sig_ht","SH");
aice = data_format(filename,"aice");
%idx = aice > 0.01;
%swh2(~idx) = 0;

w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,swh2)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('balance',31)
    %caxis([-3,3])
     title('Before')

%%

subplot(1,3,2)
filename = '/Users/noahday/Maths1/access-forcing-2010/history/iceh.2010-08.nc';



[lat,lon,row] = grid_read("gx1");
swh = data_format(filename,"wave_sig_ht"); %data_format_sector(filename,"wave_sig_ht","SH");
aice = data_format(filename,"aice");
idx = aice > 0.01;
swh(~idx) = 0;


w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
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
    title('After')


subplot(1,3,3)

w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,swh-swh2)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    cmocean('balance',31)
    caxis([-3,3])
    c.Label.String = "Signficant wave height [m]";
    title('Difference (After - Before)')





figure
subplot(1,2,1)
idx = swh>eps;
swh(~idx) = NaN;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
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
    title('New thresholds')


subplot(1,2,2)
idx = swh2>eps;
swh2(~idx) = NaN;
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,swh2)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    cmocean('balance',31)
    caxis([-3,3])
    c.Label.String = "Signficant wave height [m]";
    title('Old thresholds')


    %%

close all
    filename = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/iceh.2010-08.nc';

[lat,lon,row] = grid_read("om2");
afsdn = data_format(filename,"afsdn"); %data_format_sector(filename,"wave_sig_ht","SH");
afsd = data_format(filename,"afsd");
aice = data_format(filename,"aice");
aicen = data_format(filename,"aicen");
fsdrad = data_format(filename,"fsdrad");
%pancake = afsdn(:,:,1,1);
pancake = afsd(:,:,2);
%pancake = aicen(:,:,1);
pancake = pancake./aice;

%idx = aice > 0.01;
%swh(~idx) = 0;

figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,fsdrad)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    c = colorbar;
    cmocean('thermal')
    %caxis([0,0.3])
    %c.Label.String = 'Ice in thickest category (% of aice)';


%%
close all
 filename = '/Users/noahday/Maths1/access-forcing-2010-fixed-ic/history/iceh.2010-08.nc';
%filename = '/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh_inst.2005-01-02-00000.nc';
%filename = '/Users/noahday/GitHub/cice-dirs/runs/om2-1deg/history/iceh_ic.2005-01-01-03600.nc';
var = "fsdrad";
sector = "SH";
[lat,lon,row] = grid_read("om2");
access = data_format(filename,var);
tmask =  data_format(filename,'blkmask');
%idx = access > 0.1;
%access(~idx)= NaN;
figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
  %  setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,access)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('therm',15)
% caxis([0,10])

%
% filename = '/Users/noahday/Maths1/access-forcing-2010-no-ocean/history/iceh.2010-01.nc';
% var = "sst";
% sector = "SH";
% [lat,lon,row] = grid_read("om2");
% access = data_format(filename,var);
% tmask =  data_format(filename,'blkmask');
% %idx = access > 0.1;
% %access(~idx)= NaN;
% figure
% w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [-90 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'off')
%     setm(w, 'parallellabel', 'off')
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(lat,lon,access)
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     colorbar
%     cmocean('therm')
% %    caxis([-0,0.8])


% access = data_format(filename,var);
% fsdrad =  data_format(filename,'fsdrad');
% %idx = access > 0.1;
% %access(~idx)= NaN;
% figure
% w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [-90 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'off')
%     setm(w, 'parallellabel', 'off')
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(lat,lon,access.*fsdrad)
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     colorbar
%     cmocean('therm',15)
%  %caxis([-10,0.8])


   filename = '/Users/noahday/Maths1/spectrum_fixed/history/iceh.2007-08.nc';
   [lat,lon,row] = grid_read("gx1");
aice = data_format(filename,var);

figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    %setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,aice)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('therm',31)
%     caxis([0,1])




%%

    filename = '/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/ocean/output2010.nc';

        lat = ncread(filename,"yu_ocean");
    lon = ncread(filename, "xu_ocean");
    aice = ncread(filename,var);
    aice = aice(:,:,1) - 273.15;


figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,double(aice)')
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('therm')



%%
%filename = '/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/JRA55-do/1-4-0/8XDAILY/JRA55_03hr_forcing_2005.nc';
%lat = ncread(filename,"LAT");
%lon = ncread(filename,"LON");


filename = '/Volumes/NoahDay5TB/Gadi/JRA55-do-1-5-0/atmos/3hrPt/vas/gr/v20200916/vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_201901010000-201912312100.nc';
lat = ncread(filename,"lat");
lon = ncread(filename,"lon");
data3d = ncread(filename,"vas");
data = data3d(:,:,1);
%%
close all
figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,airtemp)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    caxis([0-30,30])

%%
filename = '/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/CAWCR/MONTHLY/2005/ww3_om2_1deg_200501.nc';

lat = ncread(filename,"LAT");
lon = ncread(filename,"LON");
airtmp3d = ncread(filename,"hs");
airtemp = airtmp3d(:,:,301);



figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,airtemp)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %cmocean('ice')
    caxis([0,5])

%%
filename = '/Users/noahday/Gadi/hs200002.nc';%'/Volumes/NoahDay5TB/raw_CAWCR/ww3_om2_1deg_200002.nc'
figure
lon = ncread(filename, 'LON');
lat = ncread(filename, 'LAT');
hs = ncread(filename, 'hs');

hs = hs(:,:,end);
%
close all
figure
w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,hs)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('haline')
    caxis([0,10])
%
%clear hs lat lon


  %%
  filename = '/Users/noahday/GitHub/random-code/ocean-2d-surface_pot_temp-1-daily-mean-ym_2013_12.nc';

  lon = ncread(filename,"xt_ocean");
  lat = ncread(filename,"yt_ocean");

  data = ncread(filename,"temp") - 273.15;

close all
figure
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data(:,:,1,1)')
    %land = shaperead('landareas', 'UseGeoCoords', true);
    %geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %cmocean('balance',31)
    %caxis([-3,3])
  %%
  filename = '/Users/noahday/Gadi/husshour1.nc';
ncdisp(filename)
  %lon = ncread(filename,"ncl1");
  %lat = ncread(filename,"ncl2");

  data = ncread(filename,"huss");% - 273.15;

 close all
 figure
surf(data)
colorbar
%%
close all
figure
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat',lon',data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    %cmocean('balance',31)
    %caxis([0,10])
%%
filename = "/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/JRA55-do/1-4-0/8XDAILY/JRA55_03hr_forcing_2018.nc";
ncdisp(filename)
lon = ncread(filename,"LON");
lat = ncread(filename,"LAT");

data = ncread(filename,"spchmd");% - 273.15;

close all
figure
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-30]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat',lon',data(:,:,1,1)')
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar

%%

filename = "/Users/noahday/Gadi/wim/1991_jan_1.nc";
%ncdisp(filename)
lon = ncread(filename,"LON");
lat = ncread(filename,"LAT");
%%
filename = "/Users/noahday/Gadi/pythoninterp/cice_t_025.nc";
ncdisp(filename)

vars = [ "airtmp" "ttlpcp" "glbrad" "dlwsfc" "spchmd" "wndewd" "wndnwd"];
filename = "/Users/noahday/Gadi/wim/2000_hour1.nc";
ncdisp(filename)
%lon = ncread(filename,"LON");
%lat = ncread(filename,"LAT");
for var = vars
    data = ncread(filename,var);% - 273.15;
    
    close all
    f = figure;
      w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-30]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'off')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 30);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        colorbar
        title(var)
   % exportgraphics(f,strcat(var,'.png'),'ContentType','image')

end

%%
filename = "/Users/noahday/Gadi/wim/iceh_inst.2002-01-01-03600.nc";
ncdisp(filename)
lon = ncread(filename,"TLON");
lat = ncread(filename,"TLAT");
%%
%clear all
filename = "/Users/noahday/Gadi/wim/iceh_inst.2003-01-01-03600.nc";%"/Users/noahday/Gadi/wim/iceh_inst.2002-01-01-03600.nc";
ncdisp(filename)
lon = ncread(filename,"TLON");
lat = ncread(filename,"TLAT");
%  Global i and j:        1134         106
var = "aice";
data = ncread(filename,var);% - 273.15;
data(1134:1134,106:106) = 10;
%idx = data < 0.01;
%data(idx) = NaN;
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-65]);
    setm(w, 'maplonlimit', [-90,90]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'on')
    setm(w, 'mlabellocation', 5);
    setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,sum(data,3))
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
    cmocean('thermal')
   % clim([-20,20])
    title(var) % "Sig. wave height [m]")%
exportgraphics(f,strcat(var,'.png'),'ContentType','image')

%%
clear all
filename = "/Users/noahday/Gadi/pythoninterp/day1.nc";
ncdisp(filename)
lon = ncread(filename,"LON");
lat = ncread(filename,"LAT");
lon = rad2deg(lon);
lat = rad2deg(lat);
%%
filename = "/Users/noahday/Gadi/wim/2000_hour1.nc";
var = "airtmp";
data = ncread(filename,var);% - 273.15;
data = data(:,:,1);
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,0]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data(:,:)')
    land = shaperead('landareas', 'UseGeoCoords', true);
%    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    colorbar
   % clim([-20,20])
    title(var)
%    exportgraphics(f,strcat(var,'.png'),'ContentType','image')    



%% Lateral melt
addpath functions
filename = "/Users/noahday/Maths1/iceh.2019-03-17.nc";
var = "dafsd_latg";
data = ncread(filename,var);% - 273.15;
NFSD = ncread(filename,"NFSD");
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
[len, wid,~] = size(data);
for i = 1:len
    for j = 1:wid
        data_cecilia(i,j) = squeeze(data(i,j,:))'*floe_binwidth';
    end
end
lon = ncread(filename,"TLON");
lat = ncread(filename,"TLAT");
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_cecilia)
    land = shaperead('landareas', 'UseGeoCoords', true);
%    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%cmocean('tempo')
    cb = colorbar;
    cb.Label.String = 'Changes to FSD [-]';
    cb.Label.Interpreter = 'latex';
    %clim([-0.05,0])
    title('After changes')
    exportgraphics(f,'after_latg.png','ContentType','image')    

filename = "/Volumes/NoahDay5TB/WIMonAlessandroRun/history/iceh.2019-03-17.nc";
var = "dafsd_latg";
data = ncread(filename,var);% - 273.15;
NFSD = ncread(filename,"NFSD");
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
[len, wid,~] = size(data);
for i = 1:len
    for j = 1:wid
        data_before(i,j) = squeeze(data(i,j,:))'*floe_binwidth';
    end
end
lon = ncread(filename,"TLON");
lat = ncread(filename,"TLAT");

f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_before)
    land = shaperead('landareas', 'UseGeoCoords', true);
%    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    cb = colorbar;
    %cmocean('tempo')
    cb.Label.String = 'Changes to FSD [-]';
    cb.Label.Interpreter = 'latex';
    %clim([-0.05,0])
    title('Before changes')
    exportgraphics(f,'before_latg.png','ContentType','image')    


f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data_cecilia-data_before)
    land = shaperead('landareas', 'UseGeoCoords', true);
%    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    cmocean('balance')
    cb = colorbar;
    clim([-0.05,0.05])
    title('Difference (after - before changes)')
    cb.Label.String = 'Changes to FSD [-]';
    cb.Label.Interpreter = 'latex';
    exportgraphics(f,'comp_latg.png','ContentType','image')    

filename = "/Volumes/NoahDay5TB/WIMonAlessandroRun/history/iceh.2019-03-17.nc";
var = "aice";
data = ncread(filename,var);% - 273.15;

f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data)
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    cmocean('ice')
    cb = colorbar;
    clim([0,1])
    title('')
    cb.Label.String = ['SIC'];
    cb.Label.Interpreter = 'latex';
    exportgraphics(f,'aice_cecilia.png','ContentType','image')

%% SWH
addpath functions
filename = "/Users/noahday/Gadi/wim/ww3_om2_025deg_20030101.nc";
var = "hs";
data = ncread(filename,var);% - 273.15;

lon = ncread(filename,"LON");
lat = ncread(filename,"LAT");
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'grid', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,data(:,:,1))
    land = shaperead('landareas', 'UseGeoCoords', true);

    %% Ocean
addpath functions
filename = '/Users/noahday/Gadi/nowaves/output2003.nc';
var = "sst";
data = ncread(filename,var);% - 273.15;
data = data(:,:,:) - 273.15;% - data(:,:,1,2);
ncdisp(filename)




%data(983:993,855:855+10) = 1;
%data(1345,969) = 100;
%tmask = ncread(filename, "tmask");
%tmask(1345,969) = 5;
%idx = data < 0.01'%== 0;
%data(idx) = NaN;
lon = ncread(filename,"xt_ocean");
lat = ncread(filename,"yt_ocean");

filename2 = '/Users/noahday/Gadi/wim/iceh.2003-01-01.nc';
lon = ncread(filename2,"TLON");
lat = ncread(filename2,"TLAT");


idx = lat < -50;
for i = 1:12
    temp = squeeze(data(:,:,i));
    temp2 = temp(idx);
    idxnan = isnan(temp2);
    data_vec(i) = sum(temp2(~idxnan));

end

figure
plot(data_vec)
%%

idx = lat>-30;
data(idx) = NaN;
%[lon lat] = meshgrid(lon_vec,lat_vec);
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', -[30,90]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', 5);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    %r = [100 75 40]
    %circlem(lat(1345,969),lon(1345,969),r,'facecolor','red','facealpha',0.5)
    pcolorm(lat,lon,data')
    colorbar
    cmocean('thermal')
    %clim([-0.5, 0.5])
    land = shaperead('landareas', 'UseGeoCoords', true);
    %title(strcat(filename(end-12:end-3)," ",var))
    title('sst Jan')
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
    exportgraphics(f,strcat(var,'6.png'),'ContentType','image')


    %% JRA55 data

filename = '/Users/noahday/Gadi/wim/JRA55_03hr_forcing_2003_5t.nc';
var = "airtmp";

for timestep = 1:5;
data = ncread(filename,var) - 273.15;
data = data(:,:,timestep);
ncdisp(filename)

%data(983:993,855:855+10) = 1;
%data(1345,969) = 100;
%tmask = ncread(filename, "tmask");
%tmask(1345,969) = 5;
%idx = data < 0.01'%== 0;
%data(idx) = NaN;
lon_vec = ncread(filename,"LON");
lat_vec = ncread(filename,"LAT");

idx = lat>-50;
data(idx) = NaN;

idx = data < -10;
data(idx) = NaN;

%[lon lat] = meshgrid(lon_vec,lat_vec);
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', -[50,90]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', 5);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    %r = [100 75 40]
    %circlem(lat(1345,969),lon(1345,969),r,'facecolor','red','facealpha',0.5)
    pcolorm(lat,lon,data)
    colorbar
    cmocean('thermal')
    clim([-2,15])
    land = shaperead('landareas', 'UseGeoCoords', true);
    title(strcat(sprintf('%g',timestep)," ",var))
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
    exportgraphics(f,strcat(var,sprintf('%g',timestep),'.png'),'ContentType','image')
    
end



%% CICE
addpath functions
clear data
filename = '/Volumes/NoahDay5TB/WIM_on/history/iceh.2011-09-30.nc'; %'/Users/noahday/Maths1/access-om2-1deg/iceh.2007-01-01.nc';

filename = '/Volumes/NoahDay/dy43/cice-dirs/runs/nowaves/history/iceh.2011-09-30.nc';

NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
var = "dafsd_newi";
data_raw = ncread(filename,var) ;%- 273.15;
[len, wid, ~] = size(data_raw);
%data = data_raw;
 for i = 1:len
    for j = 1:wid
        data(i,j) = sum(squeeze(data_raw(i,j,1)).*floe_binwidth(1)');
    end
end
% data = data(:,:,1);

%data = data.*floe_binwidth(1);
ncdisp(filename)

%data(983:993,855:855+10) = 1;
%data(1345,969) = 100;
%tmask = ncread(filename, "tmask");
%tmask(1345,969) = 5;
%idx = data < 0.01'%== 0;
%data(idx) = NaN;
lon = ncread(filename,"TLON");
lat = ncread(filename,"TLAT");

idx = lat>-30;
data(idx) = NaN;

%[lon lat] = meshgrid(lon_vec,lat_vec);
close all
conFigure(11)
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', -[30,90]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', 5);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    %r = [100 75 40]
    %circlem(lat(1345,969),lon(1345,969),r,'facecolor','red','facealpha',0.5)
    pcolorm(lat,lon,data)
    cb = colorbar;
    cb.Label.String = 'Change in FSD [-]';
    cmocean('balance') % haline
    maxval = max(max(abs(data)));
    clim([-maxval,maxval])
    %clim([-100,100])
    land = shaperead('landareas', 'UseGeoCoords', true);
    %title(strcat(filename(end-12:end-3)," ",var))
    %title(strcat("Change floe size cat. 1 ", var, " ",strcat(filename(end-12:end-3))))
    title(strcat(var," ",strcat(filename(end-12:end-3))))
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
    exportgraphics(f,strcat(var,'cice.png'),'ContentType','image')



    %% Compare two days





filename1 = '/Users/noahday/Gadi/wim/iceh.2003-01-01.nc';
filename2 = '/Users/noahday/Gadi/wim/iceh.2003-01-02.nc';
var = "uocn";
data1 = ncread(filename1,var) ;%- 273.15;
data1 = data1(:,:,1);
data2 = ncread(filename2,var) ;%- 273.15;
data2 = data2(:,:,1);
%ncdisp(filename)

%data(983:993,855:855+10) = 1;
%data(1345,969) = 100;
%tmask = ncread(filename, "tmask");
%tmask(1345,969) = 5;
%idx = data < 0.01'%== 0;
%data(idx) = NaN;
lon = ncread(filename1,"TLON");
lat = ncread(filename1,"TLAT");

idx = lat>-30;
data(idx) = NaN;

%[lon lat] = meshgrid(lon_vec,lat_vec);
close all
f = figure;
  w = worldmap('world');
    axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', -[30,90]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'on')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 30);
    setm(w, 'plabellocation', 5);
    setm(w, 'mlabelparallel', 5);
    setm(w, 'grid', 'on');
    setm(w, 'labelrotation', 'on')
    %r = [100 75 40]
    %circlem(lat(1345,969),lon(1345,969),r,'facecolor','red','facealpha',0.5)
    pcolorm(lat,lon,data1-data2)
    colorbar
    cmocean('balance')
    clim([-0.05,0.05])
    land = shaperead('landareas', 'UseGeoCoords', true);
    title(strcat(filename(end-12:end-3)," ",var))
    geoshow(w, land, 'FaceColor', [0.5 0.5 0.5])
    exportgraphics(f,strcat(var,'-diff-cice.png'),'ContentType','image')


