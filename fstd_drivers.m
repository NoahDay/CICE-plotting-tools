%% What drives the changes in the FSTD?
%% Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar
clc
% Parameters
sector = "SH";
case_name = 'ocnatmo';
filedir = '/Volumes/NoahDay5TB/cases/ocnatmo/history/iceh.';


%filedir = 'cases/ocnatmo/history/iceh.';
%'cases/momentum/history/iceh.'; %'/Volumes/NoahDay5TB/cases/momentum/history/iceh.2009-09-30.nc';
%[lat,lon,row] = grid_read(grid);
user = "a1724548";
% Make sector mask
%[len,wid] = size(lat);
fig_count = 0;


coords = sector_coords(sector);
% Take the location to be the centre of the sector
%struct('coords',sector_coords(sector),'centre_lon',(coords(1,1)+coords(2,1))/2,'north_lat',[],'south_lat',[]);
clear coords

initial_date.day = 1;
initial_date.month = 7;
initial_date.year = 2008; % 2005 is a spin-up year
if initial_date.month < 10
    initial_date.char = sprintf('%d-0%d-0%d', initial_date.year, initial_date.month, initial_date.day);
else
    initial_date.char = sprintf('%d-%d-0%d', initial_date.year, initial_date.month, initial_date.day);
end
filename = strcat(filedir,initial_date.char,'.nc');

pram.max_wind = 10; % m/s
pram.min_wind = 5; % m/s
pram.min_SIC = 0.15; % m/s
pram.total_datapoints = 1;
pram.quiver_color = 'red';
pram.label_color = 0.7*[0.4660 0.6740 0.1880];
pram.ice_edge_color = 0.7*[0.4660 0.6740 0.1880];
pram.colormap = 'ice';

[aice, sector_mask] = data_format_sector(filename,"aice","SH");
[nx,ny] = size(sector_mask);
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
info = ncinfo(filename);
attributes = info.Attributes;
coord_att = attributes(3); % Extract coordinate info
coord_string = coord_att.Value;
if coord_string(1:9) == 'TLON TLAT'
    coord_type = "t"; % T-grid
    lat = ncread(filename,"TLAT");
    lon = ncread(filename,"TLON");
else
    coord_type = "u"; % U-grid
    lat = ncread(filename,"ULAT");
    lon = ncread(filename,"ULON");
end
%filename = strcat(filedir,'2008-07','.nc');
%%
% file = strcat('iceh.','2008-07','.nc');
% fileday = 'iceh.2006-07-03.nc';
fileday = filename;
sector = "SH";
  dim=2;
 aice2 = data_format_sector(fileday,"aice",sector);
 afsdn2 = data_format_sector(fileday,"afsdn",sector);
 dafsd_weld2 = data_format_sector(fileday,"dafsd_weld",sector);
%frzmlt = data_format_sector(filename,"frzmlt",sector,dim);
 afsd2 = (sum(afsdn2,3) > eps).*1.0;
 temp_weld = 0;
 for i = 1:16
    temp_weld =  temp_weld + dafsd_weld2(:,:,i).*NFSD(i);
 end
%sum_dafsd_weld = (sum(dafsd_weld2,3)).*1.0;
%temp_data = (frzmlt>eps).*1.0;
% aminweld = 0.1;
% aice_temp = (aice2>aminweld).*1.0;
%figure(1)
  %[w1, a, output_data] = map_plot(temp_data,'frzmlt_m',sector,grid);
  figure(2)
  [w2, a, output_data] = map_plot(aice,'aice',sector,grid);
  figure(3)
  [w3, a, output_data] = map_plot(afsd2(:,:,2),'afsd',sector,grid);
  figure(4)
[w4, a, output_data] = map_plot(temp_weld,'dafsd_weld',sector,grid,[-1,10]);
% % Print all the ocean forcing


 %%
all_cat = 0;

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
    title(strcat("FSD across the ",sector," sector - ",initial_date.char))
    set(gcf,'Position',[1300 1000 800 400])
    xtickangle(45)
   
%% Impact of each term over time (integrate over space)
% This is to reproduce Figure 3 (a)-(e) in Roach et al. (2018)
datapoints = 10; % Number of days per month
date = initial_date.char;
ticker = 1;
for i = 1:datapoints
    close all
   % Get the file name
   ssd = 1;
    if ssd == 1
        filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',date,".nc");
    else
        filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
    end
%    % Aggregate data over sector - Noah Day 22/03/22 This seems wrong,
% instead apply afsd -> fsdrad conversion to calculate the change/metre
%     change_melt = aggregate_data(filename,"dafsd_latm",sector);
%     dafsd_perday(i,1) = mean(change_melt*floe_binwidth'); % Units m^2/day
%     change_growth = aggregate_data(filename,"dafsd_latg",sector);
%     dafsd_perday(i,2) = mean(change_growth*floe_binwidth'); % Units m^2/day
%     change_newi = aggregate_data(filename,"dafsd_newi",sector);
%     dafsd_perday(i,3) = mean(change_newi*floe_binwidth'); % Units m^2/day
%     change_weld = aggregate_data(filename,"dafsd_weld",sector);
%     dafsd_perday(i,4) = mean(change_weld*floe_binwidth'); % Units m^2/day
%     change_wave = aggregate_data(filename,"dafsd_wave",sector);
%     dafsd_perday(i,5) = mean(change_wave*floe_binwidth'); % Units m^2/day
% 
%    
    %change_melt = data_format_sector(filename,"dafsd_weld",sector,3);
    
    [aice, sector_mask] = data_format_sector(filename,"aice","SH");
    % Change in FSD per day (m/day)
    dafsd.latm(:,:,i) = fsd_converter(filename,"dafsd_latm","fsdrad"); 
    dafsd.latg(:,:,i) = fsd_converter(filename,"dafsd_latg","fsdrad"); 
    dafsd.newi(:,:,i) = fsd_converter(filename,"dafsd_newi","fsdrad"); 
    dafsd.weld(:,:,i) = fsd_converter(filename,"dafsd_weld","fsdrad"); 
    dafsd.wave(:,:,i) = fsd_converter(filename,"dafsd_wave","fsdrad"); 
    
    
    dafsd_SH.latm(:,:,i) = dafsd.latm(:,:,i)./sector_mask;
    dafsd_SH.latg(:,:,i) = dafsd.latg(:,:,i)./sector_mask;
    dafsd_SH.newi(:,:,i) = dafsd.newi(:,:,i)./sector_mask;
    dafsd_SH.weld(:,:,i) = dafsd.weld(:,:,i)./sector_mask;
    dafsd_SH.wave(:,:,i) = dafsd.wave(:,:,i)./sector_mask;
    
    temp.latm(:,:) = dafsd_SH.latm(:,:,i);
    temp.latg(:,:) = dafsd_SH.latg(:,:,i);
    temp.newi(:,:) = dafsd_SH.newi(:,:,i);
    temp.weld(:,:) = dafsd_SH.weld(:,:,i);
    temp.wave(:,:) = dafsd_SH.wave(:,:,i);
    
    dafsd_SH_ave.latm(i) = mean(temp.latm(~isnan(dafsd_SH.latm(:,:,i))));
    dafsd_SH_ave.latg(i) = mean(temp.latg(~isnan(dafsd_SH.latg(:,:,i))));
    dafsd_SH_ave.newi(i) = mean(temp.newi(~isnan(dafsd_SH.newi(:,:,i))));
    dafsd_SH_ave.weld(i) = mean(temp.weld(~isnan(dafsd_SH.weld(:,:,i))));
    dafsd_SH_ave.wave(i) = mean(temp.wave(~isnan(dafsd_SH.wave(:,:,i))));
    
    % Average across the sector

    
    
   % Update date
   %for j = 1:timestep
    date = update_date(date);
   %end
   if mod(i,floor(datapoints/10)) == 0
       clc
       fprintf('%g0%% complete\n',ticker);
       ticker = ticker + 1;
   end
end    

%% Plot time series
clear data_mat
close 

end_date = date;
ts.init_date = datetime(initial_date.year,initial_date.month,initial_date.day)-1;
ts.end_date = datetime(str2num(end_date(1:4)),str2num(end_date(6:7)),str2num(end_date(9:10)))-1;
ts.dates = datevec(ts.init_date:ts.end_date);
data_mat(1,:) = dafsd_SH_ave.latm; 
data_mat(2,:) = dafsd_SH_ave.latg;
data_mat(3,:) = dafsd_SH_ave.newi/1000;
data_mat(4,:) = dafsd_SH_ave.weld;
data_mat(5,:) = dafsd_SH_ave.wave/1000;
ts_wave = timeseries_plot(data_mat,strcat("Change on AFSD waves across the ",sector," sector"),'days',char(ts.init_date));



plot(ts_wave,'-', 'LineWidth',2)
    set(gcf,'Position',[1200 1000 500 200])
    ylabel('Change in L(r,h)dr')
    legend({'Lateral melt','Lateral growth','New ice','Welding','Wave induced ice fracture'},'Location','northeast')
    %ylim([-4,4])
    grid on
    xtickangle(45)
    
% Make averaged plot (monthly)    

% Get these plots compatible with Jack's plotting code

%% Plot the map over a month
% Reproduce Figure 3(k)-(o)
% sector = "SH";
% initial_date.day = 1;
% initial_date.month = 7;
% initial_date.year = 2008; % 2005 is a spin-up year
% if initial_date.month < 10
%     initial_date.char = sprintf('%d-0%d-0%d', initial_date.year, initial_date.month, initial_date.day);
% else
%     initial_date.char = sprintf('%d-%d-0%d', initial_date.year, initial_date.month, initial_date.day);
% end
% filename = strcat(filedir,initial_date.char,'.nc');
% 
% datapoints = 31; % Number of days per month
% date = initial_date.char;
% ticker = 1;
% for i = 1:datapoints
%     close all
%    % Get the file name
%    ssd = 1;
%     if ssd == 1
%         filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',date,".nc");
%     else
%         filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
%     end
%    % Aggregate data over sector
%     d_melt = data_format_sector(filename,"dafsd_latm",sector,3);
%     work = 0;
%     for j = 1:Nf 
%         work = work + d_melt(:,:,j)*floe_binwidth(j); % Units m^2/day
%     end
%     change_melt_ave(:,:,i) = work;
%     
%     d_growth = data_format_sector(filename,"dafsd_latg",sector,3);
%     work = 0;
%     for j = 1:Nf 
%         work = work + d_growth(:,:,j)*floe_binwidth(j); % Units m^2/day 
%     end
%     change_growth_ave(:,:,i) = work;
%     
%     d_newi = data_format_sector(filename,"dafsd_newi",sector,3);
%     work = 0;
%     for j = 1:Nf 
%         work = work + d_newi(:,:,j)*floe_binwidth(j); % Units m^2/day
%     end
%     change_newi_ave(:,:,i) = work;
%     
%     d_weld = data_format_sector(filename,"dafsd_weld",sector,3);
%     work = 0;
%     for j = 1:Nf 
%         work = work + d_weld(:,:,j)*floe_binwidth(j); % Units m^2/day
%     end
%     change_weld_ave(:,:,i) = work;
%     
%     d_wave = data_format_sector(filename,"dafsd_wave",sector,3);
%     work = 0;
%     for j = 1:Nf 
%         work = work + d_wave(:,:,j)*floe_binwidth(j); % Units m^2/day
%     end
%     change_wave_ave(:,:,i) = work;
% 
%     aice_daily(:,:,i) = data_format_sector(filename,"aice",sector,2);
%    % Update date
%    %for j = 1:timestep
%     date = update_date(date);
%    %end
%    if mod(i,floor(datapoints/10)) == 0
%        clc
%        fprintf('%g0%% complete\n',ticker);
%        ticker = ticker + 1;
%    end
% end    
% 
% change_melt_ave_time = sum(change_melt_ave,3)/datapoints;
% change_growth_ave_time = sum(change_growth_ave,3)/datapoints;
% change_newi_ave_time = sum(change_newi_ave,3)/datapoints;
% change_weld_ave_time = sum(change_weld_ave,3)/datapoints;
% change_wave_ave_time = sum(change_wave_ave,3)/datapoints;
% mean_aice = sum(aice_daily,3)/datapoints;

%% Plot world maps
day = 9;

plot_date = strcat(initial_date.char(1:9),sprintf('%d',day));
filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',plot_date,".nc");
close
figcount = 1;
[aice, sector_mask] = data_format_sector(filename,"aice","SH");
SIC = 0.15;
[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice,SIC,"SH",lat,lon);
%lon_ice_edge(33) = 30
[cell_lat,cell_lon] = lat_lon_finder(lat_ice_edge(33),lon_ice_edge(33),lat,lon);
line_width = 1;
tiled = tiledlayout(3,2);
tiled.Title.String = plot_date;


%tiled.YLabel.Position(1) = 2040; % change horizontal position of ylabel
%tiled.YLabel.Position(2) = 0; % change vertical position of ylabel

bounds = 15;
manual_bounds = 15;

nexttile
[figs.lat_m,a] = map_plot(dafsd_SH.latm(:,:,day),"dafsd_latm",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD lat melt"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis(1.5*[-bounds,bounds])
    

nexttile
[figs.lat_g,a] = map_plot(dafsd_SH.latg(:,:,day),"dafsd_latg",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD lat growth"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis(10^(-3)*[-bounds,bounds])
    
nexttile    
[figs.newi,a] = map_plot(dafsd_SH.newi(:,:,day),"dafsd_newi",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD new ice"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])


nexttile 
[figs.weld,a] = map_plot(dafsd_SH.weld(:,:,day),"dafsd_weld",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD welding"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-bounds,bounds])
    
% max(max(abs(change_wave_ave_time)));
nexttile
[figs.wave,a] = map_plot(dafsd_SH.wave(:,:,day),"dafsd_wave",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD wave"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    caxis([-manual_bounds,manual_bounds])
    
%% Validating dafsd

% Pick a cell
[lat_out,lon_out] = lat_lon_finder(-60,330,lat,lon);

% Take the daily difference in afsd for that cell

day1 = strcat(initial_date.char(1:9),sprintf('%d',1));
filename1 = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',day1,".nc");

day2 = strcat(initial_date.char(1:9),sprintf('%d',2));
filename2 = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',day2,".nc");


% Compare with the sum of the dafsd's

    % Get data
    afsd1 = data_format_sector(filename1,"afsd",sector,3);
    afsd2 = data_format_sector(filename2,"afsd",sector,3);
    afsd = afsd2-afsd1;
    
    latm = data_format_sector(filename2,"dafsd_latm",sector,3);
    latg = data_format_sector(filename2,"dafsd_latg",sector,3);
    newi = data_format_sector(filename2,"dafsd_newi",sector,3);
    weld = data_format_sector(filename2,"dafsd_weld",sector,3);
    wave = data_format_sector(filename2,"dafsd_wave",sector,3);
    %fsdrad = fsd_converter(filename,"afsd","fsdrad",afsd);
% Plot them
tarea =  data_format_sector(filename1,"tarea",sector,2);
tarea_cell = tarea(lon_out,lat_out);
for i = 1:Nf
    afsd_cell(i) = afsd(lon_out,lat_out,i).*tarea_cell;
    latm_cell(i) = latm(lon_out,lat_out,i).*tarea_cell;
    latg_cell(i) = latg(lon_out,lat_out,i).*tarea_cell;
    newi_cell(i) = newi(lon_out,lat_out,i).*tarea_cell;
    weld_cell(i) = weld(lon_out,lat_out,i).*tarea_cell;
    wave_cell(i) = wave(lon_out,lat_out,i).*tarea_cell;
end





figure(1)
    bar(afsd_cell)
    title('AFSD2 - AFSD1')
    
    figure(2)
    combined = [latm_cell;latg_cell;newi_cell;weld_cell;wave_cell];
    bar(combined,'stacked')
    
    
    
% [figs.wave,a] = map_plot(fsdrad,"afsd",sector);  
%     plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
%     %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
%     %t.Color = pram.label_color;
%     %set(gcf,'Position',[1200 1000 300 400])
%     a.Label.String = 'dr_a/dt (m/day)';
%     title(strcat("AFSD_2 - AFSD_1"), 'interpreter','latex','FontSize', 14)
%     cmocean('-balance',15)
%     caxis([-100,100])
%%
figcount = figcount + 1;
figure(figcount)
figs.wave = map_plot(mean_aice,"aice",sector);  
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