clear all
close all
addpath functions
%cd .. % Move up one directory
day = 6;
month = 1;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end
case_name = 'ocntest';
ticker = 1;
SIC = 0.15; 
filename = "cases/wimon/history/iceh.2005-09-01.nc";%/Volumes/NoahDay5TB/cases/forcingnowaves/history/iceh.2005-07.nc";%"prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_196201010130-196212312230.nc";
%strcat('cases/',case_name,"/history/iceh.",date,".nc");
% '/Users/noahday/Gadi/2010/JRA55_03hr_forcing_2010.nc';%'grid/gridded_ww3.glob_24m.200501.nc'; 
%filename = 'DATA/CESM/MONTHLY/ocean_forcing_clim_2D_gx1.20210330.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';


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
f3 = figure
swhdata = data_format(filename,"wave_sig_ht");
[p,a] = map_plot(swhdata,"wave_sig_ht","SH")
a.Label.String = "SWH, $m$";
a.Label.Interpreter = "latex";
exportgraphics(f3,'swh.pdf','ContentType','vector')
colormap parula

idx = swhdata > eps;
f = figure;
hist(swhdata(idx))
xlabel('Significant wave height, $m$','Interpreter','latex')
ylabel('Count','Interpreter','latex')
%exportgraphics(f,'histswh.pdf','ContentType','vector')


%% 
nbins = 50;
f = figure;
hist3([swhdata(idx), ppddata(idx)],'Nbins',[nbins,nbins],'CdataMode','auto')
xlabel('Significant wave height, $H_s$ (m)','Interpreter','latex')
    xlim([0,max(swhdata(idx))])
    ylim([0,max(ppddata(idx))])
    ylabel('Peak period, $T_p$ (s)','Interpreter','latex')
    zlabel('$r_{max}$','Interpreter','latex')
    s.EdgeColor = 'none';
    a = colorbar;
    a.Label.String = 'Counts';
    a.Label.Interpreter = 'latex';
    a.TickLabelInterpreter = 'latex';
    C=jet(15);
    
    colormap(C)
    view(0,90)
    
exportgraphics(f,'countswhppd.pdf','ContentType','image')
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
