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
filename = "cases/31freq/history/iceh.2008-07-05.nc";%/Volumes/NoahDay5TB/cases/forcingnowaves/history/iceh.2005-07.nc";%"prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_196201010130-196212312230.nc";
%strcat('cases/',case_name,"/history/iceh.",date,".nc");
% '/Users/noahday/Gadi/2010/JRA55_03hr_forcing_2010.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
%filename = 'DATA/CESM/MONTHLY/ocean_forcing_clim_2D_gx1.20210330.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';

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
f = figure;
dafsd3d = data_format(filename,"dafsd_newi");
t = tiledlayout(4,4);
t.TileSpacing = 'compact';
for i = 1:16
    nexttile
    dafsd = dafsd3d(:,:,i);
    idx = dafsd < eps;
    dafsd(idx) = NaN;
    if isnan(max(max(dafsd)))
        maxi = 1;
    else
        maxi = max(max(dafsd));
    end
    [p,a] = map_plot(dafsd,"dafsd_newi","SH","gx1",[0,maxi]);
    if isnan(max(max(dafsd)))
        title('No Change')
    else
        maxi = max(max(dafsd));
    end
    %a.Label.String = "Change in FSD";
    a.Label.Interpreter = "latex";
   % exportgraphics(f,'dafsd.pdf','ContentType','vector')
    colormap parula
end
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
