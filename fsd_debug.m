% Debug the evolution of the FSD
%% 1. Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar
clc
% Parameters
sector = "SH";
grid = 'gx1';
filedir = '/Users/noahday/Maths1/fsd_debug/withoutWIM/iceh_inst.';
%filedir = 'cases/ocnatmo/history/iceh.';
%'cases/momentum/history/iceh.'; %'/Volumes/NoahDay5TB/cases/momentum/history/iceh.2009-09-30.nc';
[lat,lon,row] = grid_read(grid);
user = "a1724548";
% Make sector mask
%[len,wid] = size(lat);
fig_count = 0;


coords = sector_coords(sector);
% Take the location to be the centre of the sector
%struct('coords',sector_coords(sector),'centre_lon',(coords(1,1)+coords(2,1))/2,'north_lat',[],'south_lat',[]);
clear coords
initial_date.hour = 3600;
for i = 1:23
    initial_date.day = 1;
    initial_date.month = 7;
    initial_date.year = 2008; % 2005 is a spin-up year
    if initial_date.hour < 10000
        initial_date.char = sprintf('%d-0%d-0%d-0%d', initial_date.year, initial_date.month, initial_date.day,initial_date.hour);
    else
        initial_date.char = sprintf('%d-0%d-0%d-%d', initial_date.year, initial_date.month, initial_date.day,initial_date.hour);
    end
    filenames(i,:) = strcat(filedir,initial_date.char,'.nc');
    initial_date.hour = initial_date.hour + 3600;
end

pram.max_wind = 10; % m/s
pram.min_wind = 5; % m/s
pram.min_SIC = 0.15; % m/s
pram.total_datapoints = 1;
pram.quiver_color = 'red';
pram.label_color = 0.7*[0.4660 0.6740 0.1880];
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
pram.colormap = 'ice';

file_ic = '/Users/noahday/Maths1/fsd_debug/withWIM/iceh_ic.2008-07-01-03600.nc';
filenameday = '/Users/noahday/Maths1/fsd_debug/withWIM/iceh.2008-07-01.nc';
[aice, sector_mask] = data_format_sector(filenames(i,:),"aice",sector);
[nx,ny] = size(sector_mask);
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

NFSD = ncread(filenameday,"NFSD");
NCAT = ncread(filenameday,"NCAT");
%filename = strcat(filedir,'2008-07','.nc');
 %%
all_cat = 0;
filename = filenames(1,:);
%  Map by term for all categories
if all_cat == 1
    close all
    fig_count = 0;
    for i = 1:16
        fig_count = fig_count+1;
        figure(fig_count)
        plot_map_fsd(filename,"dafsd_wave", i, sector)
    end
end

% FSD of the SA sector
fsd = fsd_converter(filename,"afsdn","fsd");
count = 1;
for i = 1:nx
    for j = 1:ny
        if aice(i,j) > eps % Only take the cells where there IS ice
            for n = 1:Nf
                fsd_tab(count,n) = fsd(i,j,n)./aice(i,j);
            end
            aice_tab(count) = aice(i,j);
            count = count + 1;
        end
    end
end

%
cum_fsd = sum(fsd_tab);
normalised_fsd = cum_fsd./sum(cum_fsd);
% FSD histogram
fig_count = fig_count + 1;
figure(fig_count)
bar(1:Nf,normalised_fsd,'hist')
xlabel('FSD radius (m)')
ylabel('Fraction of sea ice area (normalised)')
Nf = numel(NFSD);
format shortg
r_NFSD = round(round(floe_rad_h,1));
r1_NFSD = round(round(floe_rad_l,1));
s_NFSD = num2str(r_NFSD);
for i = 1:Nf
    lab{i} = sprintf('[%g, %g]',r1_NFSD(i),r_NFSD(i));
end
    xticks(1:Nf)
    xticklabels(lab)
    title(strcat("FSD across the ",sector, ", ",filename(end-18:end-3)," without waves"))
    set(gcf,'Position',[1300 1000 800 400])
    xtickangle(45)
   
%% 2. Impact of each term over time (integrate over space)
% This is to reproduce Figure 3 (a)-(e) in Roach et al. (2018)
datapoints = 23; % Number of days per month
%date = initial_date.char;
ticker = 1;
for i = 1:datapoints
    close all
   % Get the file name
%    ssd = 1;
%     if ssd == 1
%         filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',date,".nc");
%     else
%         filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
%     end

   filename = filenames(i,:);
   % Aggregate data over sector
    change_melt = aggregate_data(filename,"dafsd_latm",sector);
    dafsd_perday(i,1) = mean(change_melt*floe_binwidth'); % Units m^2/day
    change_growth = aggregate_data(filename,"dafsd_latg",sector);
    dafsd_perday(i,2) = mean(change_growth*floe_binwidth'); % Units m^2/day
    change_newi = aggregate_data(filename,"dafsd_newi",sector);
    dafsd_perday(i,3) = mean(change_newi*floe_binwidth'); % Units m^2/day
    change_weld = aggregate_data(filename,"dafsd_weld",sector);
    dafsd_perday(i,4) = mean(change_weld*floe_binwidth'); % Units m^2/day
    change_wave = aggregate_data(filename,"dafsd_wave",sector);
    dafsd_perday(i,5) = mean(change_wave*floe_binwidth'); % Units m^2/day

   
   
   % Update date
   %for j = 1:timestep
    %date = update_date(date);
   %end

   if mod(i,floor(datapoints/10)) == 0
       clc
       fprintf('%g0%% complete\n',ticker);
       ticker = ticker + 1;
   end
end    

% Plot time series
end_date = date;
ts.init_date = datetime(initial_date.year,initial_date.month,initial_date.day)-1;
ts.end_date = datetime(str2num(end_date(1:4)),str2num(end_date(6:7)),str2num(end_date(9:10)))-1;
ts.dates = datevec(ts.init_date:ts.end_date);
ts_wave = timeseries_plot(dafsd_perday,strcat("Change on AFSD waves across the ",sector," sector"),'days',char(ts.init_date));



plot(ts_wave,'-', 'LineWidth',2)
    set(gcf,'Position',[1200 1000 500 200])
    ylabel('Change in L(r,h)dr')
    legend({'Lateral melt','Lateral growth','New ice','Welding','Wave induced ice fracture'},'Location','northeast')
    %ylim([-4,4])
    grid on
    xtickangle(45)
    
% Make averaged plot (monthly)    

% Get these plots compatible with Jack's plotting code

%% Plot the maps
% Reproduce Figure 3(k)-(o)
clear filename d_melt d_growth d_newi d_weld d_wave change_wave change_weld change_newi change_growth change_melt
clear change_melt_ave_time change_growth_ave_time change_newi_ave_time change_weld_ave_time change_wave_ave_time
clear change_growth_ave change_melt_ave change_newi_ave change_weld_ave change_wave_ave
sector = "SH";

%filename = strcat(filedir,initial_date.char,'.nc');

datapoints = 1; % Number of hours in a day
%date = initial_date.char;
ticker = 1;


for i = 1%1:datapoints
    close all
   % Get the file name
   %filename = filenames(23,:);
   filename = ['/Users/noahday/Maths1/fsd_debug/withoutWIM/iceh_inst.2008-07-03-00000.nc'];
   % Aggregate data over sector
    d_melt = data_format_sector(filename,"dafsd_latm",sector,3);
    work = 0;
    for j = 1:Nf 
        work = work + d_melt(:,:,j)*floe_binwidth(j); % Units m^2/day
    end
    change_melt_ave(:,:,i) = work;
    clear work
    d_growth = data_format_sector(filename,"dafsd_latg",sector,3);
    work = 0;
    for j = 1:Nf 
        work = work + d_growth(:,:,j)*floe_binwidth(j); % Units m^2/day 
    end
    change_growth_ave(:,:,i) = work;
    clear work
    d_newi = data_format_sector(filename,"dafsd_newi",sector,3);
    work = 0;
    for j = 1:Nf 
        work = work + d_newi(:,:,j)*floe_binwidth(j); % Units m^2/day
    end
    change_newi_ave(:,:,i) = work;
    clear work

    d_weld = data_format_sector(filename,"dafsd_weld",sector,3);
    work = 0;
    for j = 1:Nf 
        work = work + d_weld(:,:,j)*floe_binwidth(j); % Units m^2/day
    end
    change_weld_ave(:,:,i) = work;
    clear work

    d_wave = data_format_sector(filename,"dafsd_wave",sector,3);
    work = 0;
    for j = 1:Nf 
        work = work + d_wave(:,:,j)*floe_binwidth(j); % Units m^2/day
    end
    change_wave_ave(:,:,i) = work;
    clear work
    aice_daily(:,:,i) = data_format_sector(filename,"aice",sector,2);
   % Update date
   %for j = 1:timestep
    %date = update_date(date);
   %end
   if mod(i,floor(datapoints/10)) == 0
       clc
       fprintf('%g0%% complete\n',ticker);
       ticker = ticker + 1;
   end
end    

change_melt_ave_time = sum(change_melt_ave,3)/datapoints;
change_growth_ave_time = sum(change_growth_ave,3)/datapoints;
change_newi_ave_time = sum(change_newi_ave,3)/datapoints;
change_weld_ave_time = sum(change_weld_ave,3)/datapoints;
change_wave_ave_time = sum(change_wave_ave,3)/datapoints;
mean_aice = sum(aice_daily,3)/datapoints;

SWH = data_format_sector(filename,"wave_sig_ht",sector,2);
%for i = 1:datapoints
%    temp = change_weld_ave(:,:,i);
%   ts(i) =  sum(temp(~isnan(temp)));
%end
% Plot world maps
close all
figcount = 1;
[aice, sector_mask] = data_format_sector(filename,"aice","SH");
SIC = 0.15;
[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,"SH",lat,lon);
%lon_ice_edge(33) = 30 
[cell_lat,cell_lon] = lat_lon_finder(lat_ice_edge(33),lon_ice_edge(33),lat,lon);
line_width = 1;
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
plot_time = filename(end-18:end-3);
t.Title.String = plot_time;
t.Title.Interpreter = 'latex';
 %set(tiled,'Position',[1200 1000 300 400])
bounds = max(max(abs(change_melt_ave_time)));
manual_bounds = 15;


nexttile
figs.lat_m = map_plot(change_melt_ave_time,"dafsd_latm",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Change in FSD lat melt"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])
    
bounds = max(max(abs(change_growth_ave_time)));
nexttile
figs.lat_g = map_plot(change_growth_ave_time,"dafsd_latg",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Change in FSD lat growth"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])

bounds = max(max(abs(change_newi_ave_time)));    
nexttile    
figs.newi = map_plot(change_newi_ave_time,"dafsd_newi",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Change in FSD new ice"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])

bounds = max(max(abs(change_weld_ave_time)));
nexttile 
figs.weld = map_plot(change_weld_ave_time,"dafsd_weld",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Change in FSD welding"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])
    
% max(max(abs(change_wave_ave_time)));
nexttile
figs.wave = map_plot(change_wave_ave_time,"dafsd_wave",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Change in FSD wave"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])

nexttile
figs.swh = map_plot(SWH,"wave_sig_ht",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    title(strcat("SWH"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])

%%
figcount = figcount + 1;
figure(figcount)
figs.wave = map_plot(mean_aice,"aice",sector);  

%% Evolution of the FSD over a single timestep
file_ic = '/Users/noahday/Maths1/fsd_debug/withWIM/iceh_ic.2008-07-01-03600.nc';
file_wim = '/Users/noahday/Maths1/fsd_debug/withWIM/iceh.2008-07-01-03600.nc';
file_nowim = '/Users/noahday/Maths1/fsd_debug/withoutWIM/iceh.2008-07-01-03600.nc';

% FSD of the SA sector
fsd = fsd_converter(filename,"afsdn","fsd");
count = 1;
for i = 1:nx
    for j = 1:ny
        if aice(i,j) > eps % Only take the cells where there IS ice
            for n = 1:Nf
                fsd_tab(count,n) = fsd(i,j,n)./aice(i,j);
            end
            aice_tab(count) = aice(i,j);
            count = count + 1;
        end
    end
end

%
cum_fsd = sum(fsd_tab);
normalised_fsd = cum_fsd./sum(cum_fsd);
% FSD histogram
fig_count = fig_count + 1;
figure(fig_count)
bar(1:Nf,normalised_fsd,'hist')
xlabel('FSD radius (m)')
ylabel('Fraction of sea ice area (normalised)')
Nf = numel(NFSD);
format shortg
r_NFSD = round(round(floe_rad_h,1));
r1_NFSD = round(round(floe_rad_l,1));
s_NFSD = num2str(r_NFSD);
for i = 1:Nf
    lab{i} = sprintf('[%g, %g]',r1_NFSD(i),r_NFSD(i));
end
    xticks(1:Nf)
    xticklabels(lab)
    title(strcat("FSD across the ",sector, ", ",filename(end-18:end-3)," without waves"))
    set(gcf,'Position',[1300 1000 800 400])
    xtickangle(45)
%% Validate the FSD converter   
sector = "SH";
fsdrad_test = fsd_converter(filename,"afsdn","fsdrad");
fsdrad_actual = data_format_sector(filename,"fsdrad",sector);

[lat_out,lon_out] = lat_lon_finder(-65,30,lat,lon);


fsdrad_test(lon_out,lat_out)
fsdrad_actual(lon_out,lat_out)

dafsd_test = fsd_converter(filename,"dafsd_wel","fsdrad");
dafsd_test(lon_out,lat_out)/1000;

%% Functions
function map = plot_map_fsd(filename,variable,cat,sector)
    % Define ice edge 
    [lat,lon,~] = grid_read("gx1");
    aice_data = data_format_sector(filename,"aice",sector);
    idx = aice_data < eps;
    aice_data(idx) = NaN;
    lat_vec = reshape(lat,1,[]);
    lon_vec = reshape(lon,1,[]);

    SIC = 0.15;
    [lat_ice_edge, lon_ice_edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);
    
    % Get fsd data
    data3D = data_format_sector(filename,variable,sector,3);
    data = data3D(:,:,cat);
    idx = data < eps;
    data(idx) = NaN;
    fig = map_plot(data,variable,sector);  
        plotm(lat_ice_edge,lon_ice_edge,'-','color','k','LineWidth',2)
        t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',SIC),'HorizontalAlignment','right');
        set(gcf,'Position',[1200 1000 300 400])
        title(sprintf("Change in FSD due to lat melt - cat %g",cat), 'interpreter','latex','FontSize', 14)
        colormap(parula(10))


end

function output = aggregate_data(filename,variable,sector)
    data = data_format_sector(filename,variable,sector,3);
    aice = data_format_sector(filename,"aice",sector);
    NFSD = ncread(filename,"NFSD");
    Nf = numel(NFSD);
    [nx, ny] = size(aice);
    count = 1;
    for i = 1:nx
        for j = 1:ny
            if aice(i,j) > eps % Only take the cells where there IS ice
                for n = 1:Nf
                    dafsd_tab(count,n) = data(i,j,n);
                end
                %aice_tab(count) = aice(i,j);
                count = count + 1;
            else
                dafsd_tab(count,:) = 0;
            end
        end
    end
    output = dafsd_tab;
end

function ts = timeseries_plot(data,name,units,startdate)
    ts = timeseries(data,1:length(data));
    ts.Name = name;
    ts.TimeInfo.Units = units;
    ts.TimeInfo.StartDate = startdate;     % Set start date.
    ts.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.
end