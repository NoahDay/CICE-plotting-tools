%% What drives the changes in the FSTD?
%% Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar
clc
% Parameters
sector = "SH";
case_name = 'wimon';
filedir = 'cases/wimon/history/iceh.';


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
initial_date.month = 9;
initial_date.year = 2005; % 2005 is a spin-up year
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


%[aice, sector_mask] = data_format_sector(filename,"aice","SH");
nx = 321;
ny = 384;
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
info = ncinfo(filename,"afsdn");
attributes = info.Attributes;
coord_att = attributes(3); % Extract coordinate info
coord_string = coord_att.Value;
if coord_string(1:9) == 'TLON TLAT'
    coord_type = "t"; % T-grid
    lat = ncread(filename,"TLAT");
    lon = ncread(filename,"TLON");
    dim_info = ncinfo(filename,"TLAT");
else
    coord_type = "u"; % U-grid
    lat = ncread(filename,"ULAT");
    lon = ncread(filename,"ULON");
end
data_size = dim_info.Size;
dim = length(data_size);
if data_size(1) == 320 && data_size(2) == 384
    grid = "gx1";
    row = 37;
    lat = rearrange_matrix(lat,row,2);
    lon = rearrange_matrix(lon,row,2);

    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
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
  [w2, a, output_data] = map_plot(aice2,'aice',sector,grid);
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
clear temp dafsd dafsd_SH data
clc
cases = ["forcingoff","wimon","wimoff"];%["profile","nowaves"];
datapoints = 24; % Number of days per month
date = initial_date.char;
ssd = 0;
sector = "SH";
[data] = read_in_fsd_data(cases,date,datapoints,sector,ssd,"d");

%% Plot time series
clear data_mat
close 

end_date = date;
ts.init_date = datetime(initial_date.year,initial_date.month,initial_date.day)-1;
ts.end_date = datetime(str2num(end_date(1:4)),str2num(end_date(6:7)),str2num(end_date(9:10)))-1;
%ts.end_date = datetime(2005,08,01);
ts.dates = datevec(ts.init_date:ts.end_date);
data_mat(1,:) = data.latm.ra(:,1); 
data_mat(2,:) = data.latg.ra(:,1);
data_mat(3,:) = data.newi.ra(:,1);
data_mat(4,:) = data.weld.ra(:,1);
data_mat(5,:) = data.wave.ra(:,1);
ts_wave = timeseries_plot(data_mat,strcat("Wave forcing = 0"),'days',char(ts.init_date));



plot(ts_wave,'-', 'LineWidth',2)
    set(gcf,'Position',[1200 1000 600 300])
    ylabel('Change in $r_a$ (m/day)','Interpreter','Latex')
    legend({'Lateral melt','Lateral growth','New ice','Welding','Wave induced ice fracture'},'Location','northeast')
    ylim([-0.6,1.4])
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
%set(gca,'ColorScale','log')

nexttile
cust_bounds =  max(max(abs(dafsd_SH.latm(:,:,day))));
[Ticks,TickLabels] = log_ticks(cust_bounds,10,'balance');

[figs.lat_m,a] = map_plot(1/10*(dafsd_SH.latm(:,:,day)),"dafsd_latm",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD lat melt"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    a.Ticks = Ticks;
    a.TickLabels = TickLabels;
   caxis([-log(cust_bounds) log(cust_bounds)]);

nexttile
cust_bounds =  max(max(abs(dafsd_SH.latg(:,:,day))));
%[Ticks,TickLabels] = log_ticks(cust_bounds,10,'balance');

[figs.lat_g,a] = map_plot(dafsd_SH.latg(:,:,day),"dafsd_latg",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD lat growth"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    %a.Ticks = Ticks;
   % a.TickLabels = TickLabels;
   caxis([-cust_bounds cust_bounds]);
    
nexttile    
cust_bounds =  max(max(abs(dafsd_SH.newi(:,:,day))));
[Ticks,TickLabels] = log_ticks(cust_bounds,10,'balance');

[figs.newi,a] = map_plot(dafsd_SH.newi(:,:,day),"dafsd_newi",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD new ice"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    a.Ticks = Ticks;
    a.TickLabels = TickLabels;
   caxis([-log(cust_bounds) log(cust_bounds)]);


nexttile 
cust_bounds =  max(max(abs(dafsd_SH.weld(:,:,day))));
[figs.weld,a] = map_plot(dafsd_SH.weld(:,:,day),"dafsd_weld",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD welding"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    a.Ticks = Ticks;
    a.TickLabels = TickLabels;
   caxis([-log(cust_bounds) log(cust_bounds)]);
    
cust_bounds =  max(max(abs(dafsd_SH.wave(:,:,day))));
[Ticks,TickLabels] = log_ticks(cust_bounds,10,'balance');

nexttile
[figs.wave,a] = map_plot(dafsd_SH.wave(:,:,day),"dafsd_wave",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',line_width)
    %t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    %t.Color = pram.label_color;
    %set(gcf,'Position',[1200 1000 300 400])
    a.Label.String = 'dr_a/dt (m/day)';
    title(strcat("Change in FSD wave"), 'interpreter','latex','FontSize', 14)
    cmocean('-balance',15)
    a.Ticks = Ticks;
    a.TickLabels = TickLabels;
   caxis([-log(cust_bounds) log(cust_bounds)]);
    
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
    afsd1 = data_format_sector(filename1,"afsd",sector);
    afsd2 = data_format_sector(filename2,"afsd",sector);
    afsd = afsd2-afsd1;
    
    latm = data_format_sector(filename2,"dafsd_latm",sector);
    latg = data_format_sector(filename2,"dafsd_latg",sector);
    newi = data_format_sector(filename2,"dafsd_newi",sector);
    weld = data_format_sector(filename2,"dafsd_weld",sector);
    wave = data_format_sector(filename2,"dafsd_wave",sector);
    %fsdrad = fsd_converter(filename,"afsd","fsdrad",afsd);
% Plot them
tarea =  data_format_sector(filename1,"tarea",sector);
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

%% Plot the changes of each floe size category over time
%% Plotting unscaled dafsd
NFSD = ncread(filename,'NFSD');
Nf = length(NFSD);
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
f_bin_width = floe_rad_h - floe_rad_l;

r_NFSD = round(round(floe_rad_h,1));
r1_NFSD = round(round(floe_rad_l,1));
s_NFSD = num2str(r_NFSD);
for i = 1:Nf
    lab{i} = sprintf('[%g, %g]',r1_NFSD(i),r_NFSD(i));
end


f = figure(1);
% Colors: latg, latm, newi, weld, wave
colors = [[0.4660 0.6740 0.1880]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [0 0.4470 0.7410]];
t = tiledlayout(2,8);
t.TileSpacing = 'compact';
for nf = 1:16
    cat = i;
    nexttile
    dafsd_data = [data.latm.ave(nf,:,1); data.latg.ave(nf,:,1); data.newi.ave(nf,:,1); data.weld.ave(nf,:,1); data.wave.ave(nf,:,1)].*f_bin_width(nf);%    dafsd_data2 = [data.latm.ave(nf,:,2); data.latg.ave(nf,:,2); data.newi.ave(nf,:,2); data.weld.ave(nf,:,2); data.wave.ave(nf,:,2)].*f_bin_width(nf);
    dafsd_data2 = [data.latm.ave(nf,:,2); data.latg.ave(nf,:,2); data.newi.ave(nf,:,2); data.weld.ave(nf,:,2); data.wave.ave(nf,:,2)].*f_bin_width(nf);
    %dataoff = [data.off.latg(:,cat), data.off.latm(:,cat), data.off.newi(:,cat), data.off.weld(:,cat), data.off.wave(:,cat)];
    [len,wid] = size(dafsd_data);
    hold on
    for j = 1:len
        p(j) = plot(1:wid, dafsd_data(j,:), '--','LineWidth',2);
        k(j) = plot(1:wid, dafsd_data2(j,:), ':*','LineWidth',2);
        set(p(j),'Color',colors(j,:));
        set(k(j),'Color',colors(j,:));
    end
    hold off
    grid on
    if nf == 1
        legend({'lat melt','lat growth','new ice','weld', 'wave'},'Position',[0.94 0.4 0.05 0.3],'FontSize',12)
        legend boxoff 
        %ylim([-0.8,0.8])
        %title(strcat(lab(i)));%title(sprintf('FSD cat %d (different y-limits!)',cat))
        %xlabel("(different y-limits!)")
    %elseif i == 9
        %ylim([-0.3,0.3])
       % title(lab(i));%sprintf('FSD cat %d',cat))
    else
        %ylim([-0.3,0.3])
        title(lab(i));%title(sprintf('FSD cat %d',cat))
    end
    title(lab(nf))
    maximum = max(max(abs(dafsd_data(:,:))));
    if nf < 13 
        ylim([-0.1,0.1])
    elseif nf < 16
        ylim([-0.5,0.5])
    else
        ylim([-5,5])
    end%
    %xlim([1,5])
%     if isnan(data(j,:)) 
%         ylim([-0.3,0.3]);
%     elseif maximum > 1
%         ylim([-5,5]);
%     elseif maximum > 0.1
%         ylim([-0.5,0.5]);
%     elseif maximum > 0.01
%         ylim([-0.1,0.1]);
%     elseif maximum > 0.001
%         ylim([-0.005,0.005]);
%     elseif maximum > 0.0001
%         ylim([-0.0005,0.0005]);
%     elseif maximum > 0.00001
%         ylim([-0.00005,0.00005]);
%     else
%          ylim([-0.000005,0.000005]);
%     end
    %ylim([-0.001,0.001])
%     if i == 1 || i == 2 || i == 9
%     else
%         set(gca,'Yticklabel',[]) 
%         %set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
%     end
%     if sum(i == 1:8) == 1
%         set(gca,'Xticklabel',[]) 
%    end
    %clear data
end
xlabel(t,'Months','FontSize',16,'Interpreter','Latex')
ylabel(t,'df(r)dr/dt (conc/day)','FontSize',16,'Interpreter','Latex')
title(t,'Change in FSD averaged across all of Antarctica','FontSize',16,'Interpreter','Latex')
f.Position = [100 100 1200 400];

%% Plotting unscaled dafsd scaled for afsd
close 
NFSD = ncread(filename,'NFSD');
Nf = length(NFSD);
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
f_bin_width = floe_rad_h - floe_rad_l;
floe_binwidth = f_bin_width;
r_NFSD = round(round(floe_rad_h,1));
r1_NFSD = round(round(floe_rad_l,1));
s_NFSD = num2str(r_NFSD);
for i = 1:Nf
    lab{i} = sprintf('[%g, %g]',r1_NFSD(i),r_NFSD(i));
end


f = figure(1);
% Colors: latg, latm, newi, weld, wave
colors = [[0.4660 0.6740 0.1880]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [0 0.4470 0.7410]];
t = tiledlayout(2,8);
t.TileSpacing = 'compact';
for nf = 1:16
    cat = i;
    nexttile
    dafsd = [data.latm.ave(nf,:,1); data.latg.ave(nf,:,1); data.newi.ave(nf,:,1); data.weld.ave(nf,:,1); data.wave.ave(nf,:,1)]./data.afsd.ave(nf,:,1);
    %dataoff = [data.off.latg(:,cat), data.off.latm(:,cat), data.off.newi(:,cat), data.off.weld(:,cat), data.off.wave(:,cat)];
    [len,wid] = size(dafsd);
    hold on
    for j = 1:len
        p(j) = plot(1:wid, dafsd(j,:), '--','LineWidth',2);
       % k(j) = plot(1:6, dataoff(:,j), '-', 'LineWidth',2);
        set(p(j),'Color',colors(j,:));
       % set(k(j),'Color',colors(j,:));
    end
    grid on
    if nf == 1
        legend({'latm','latg','newi','weld', 'wave'},'Position',[0.94 0.4 0.05 0.3])
        %ylim([-0.8,0.8])
        %title(strcat(lab(i)));%title(sprintf('FSD cat %d (different y-limits!)',cat))
        %xlabel("(different y-limits!)")
    %elseif i == 9
        %ylim([-0.3,0.3])
       % title(lab(i));%sprintf('FSD cat %d',cat))
    else
        %ylim([-0.3,0.3])
        title(lab(i));%title(sprintf('FSD cat %d',cat))
    end
    title(lab(nf))

    maximum = max(abs(dafsd(j,:)));
    %if isnan(dafsd(j,:))
    %    ylim([-0.3,0.3]);
    %else
    %    ylim([-2*maximum,2*maximum]);
    %end
    clear maximum
    %if nf < 13
    %    ylim([-0.1,0.1])
    %elseif nf < 16
    %    ylim([-0.5,0.5])
    %else
    %    ylim([-5,5])
    %end
%     if i == 1 || i == 2 || i == 9
%     else
%         set(gca,'Yticklabel',[]) 
%         %set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
%     end
%     if sum(i == 1:8) == 1
%         set(gca,'Xticklabel',[]) 
%    end

end
xlabel(t,'time steps','FontSize',14)
ylabel(t,'dafsd/afsd per day (x/day)','FontSize',14)
title(t,'DAFSD across all of Antarctica scaled for AFSD')
f.Position = [100 100 1200 400];
hold off

%% Raw AFSD
close
% ax = gca;
% ax.Title.String = '';
% ax.Title.Interpreter = 'latex';
% set(groot,'DefaultAxesTitle', ax.Title);

f2 = figure(2);
t2 = tiledlayout(2,6);
t2.TileSpacing = 'compact';


for i = 1:datapoints
    cat = i;
    nexttile
    hold on
    bar_data = [data.afsd.ave(:,i,2)];%,data.afsd.ave(:,i,2),data.afsd.ave(:,i,2)-data.afsd.ave(:,i,1)];
    bar(bar_data);
    title(sprintf("Month %d",i))
    if i == 1
        legend('Forcing off','case2','Difference')
    end
   % Nf = numel(NFSD);
   % format shortg
   % r_NFSD = round(round(floe_rad_h,1));
   % r1_NFSD = round(round(floe_rad_l,1));
   % s_NFSD = num2str(r_NFSD);
   % for i = 1:Nf
   %     lab{i} = sprintf('[%g, %g]',r1_NFSD(i),r_NFSD(i));
   % end
   % xticks(1:Nf)
   % xticklabels(lab)
   % xtickangle(45)
   clear bar_data
end


   
xlabel(t2,'FSD cats','FontSize',14)
ylabel(t2,'afsd','FontSize',14)
title(t2,'Southern hemisphere FSD - afsd, F(r)')
f2.Position = [100 100 1200 400];
hold off

% F(r)dr
f3 = figure(3);
t3 = tiledlayout(2,6);
t3.TileSpacing = 'compact';


for i = 1:datapoints
    cat = i;
    nexttile
    %[len,wid] = size(data);
    hold on
    bar_data = [data.afsd.ave(:,i,2).*f_bin_width'];%,data.afsd.ave(:,i,2).*floe_binwidth',(data.afsd.ave(:,i,2)-data.afsd.ave(:,i,1)).*floe_binwidth'];
    bar(bar_data);
     if i == 1
        legend('Forcing off','case2','Difference')
    end
    title(sprintf("Month %d",i))
    clear bar_data
end

xlabel(t3,'FSD cats','FontSize',14)
ylabel(t3,'FSD','FontSize',14)
title(t3,'Southern hemisphere FSD, F(r)dr')
f3.Position = [100 100 1200 400];
hold off

%% AFSD converted to number FSD
f4 = figure(4);
t4 = tiledlayout(2,6);
t4.TileSpacing = 'compact';
alpha = 0.66;
%number_fsd = data_fsd./(4*alpha)

for i = 1:datapoints
    cat = i;
    nexttile
    %[len,wid] = size(data);
    hold on
    num_data = [data.afsd.ave(:,i,2).*floe_binwidth'./(4*alpha*NFSD.^2)];%, data.afsd.ave(:,i,2).*floe_binwidth'./(4*alpha*NFSD)];
    bar(num_data);
    title(sprintf("Month %d",i))
     if i == 1
        legend('Forcing on waves off','case2','Difference')
    end
    clear num_data
end

xlabel(t4,'FSD cats','FontSize',14)
ylabel(t4,'Number FSD','FontSize',14)
title4 = title(t4,'Southern hemisphere number FSD, F^N(r)dr');
%titl4.Interpreter = 'Latex'
f4.Position = [100 100 1200 400];
hold off

%% Reproducing Roach et al. (2018) Fig 2. Hemishphere FSD
close all
clc
datapoints = 12;
roach_data_sep = [0.3*10^4, 3*10^1, 10^0, 0.8*10^(-1), 0.9*10^(-2), 1.2*10^(-3), 0.5*10^(-4), 0.7*10^(-5), 10^(-6), 1.3*10^(-7), 10^(-7), 10^(-4)];
f5 = figure(5);
t5 = tiledlayout(1,1);%(5,2);
t5.TileSpacing = 'compact';
alpha = 0.66;
cust_bounds =  max(NFSD);
xtick = 10.^(0:4);
xticklab = cellstr(num2str(round(log10(xtick(:))), '10^{%d}'));
ytick = 10.^(-30:2:5);
yticklab = cellstr(num2str(round(log10(ytick(:))), '10^{%d}'));

% WIM ON
filename = 'cases/wimon/history/iceh.2005-09-30.nc';
case1 = data_format_sector(filename,"afsdn",sector);
%n_fsd = fsd_converter(filename,"afsdn","fsdrad",sector);

raw_data = case1;

alpha = 0.66; % Rothrock and Thorndike (1984)
% Get tracer array
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
afsdn = data_format_sector(filename,"afsdn",sector);
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                trcrn(i,j,k,n) = raw_data(i,j,k,n);%.*floe_binwidth(k)./aicen_data(i,j,n); % something x Length: (something m)
            end
        end
    end
end
% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*floe_rad_c.^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                nfstd(i,j,k,n) = trcrn(i,j,k,n)/floe_area_c(k); % Dimensionless: something / Length
            end
        end
    end
end
% Intergrate nfstd to get nfsd

for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            processed_data(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
        end
    end
end

icemask = aice2 > 0.01;
data1(:,:,:) = processed_data;

for i = 1:Nf
    temp = data1(:,:,i);
    data1_overspace(i) = mean(temp(icemask),'all');
end


% WIM off
filename = 'cases/wimoff/history/iceh.2005-09-30.nc';
case1 = data_format_sector(filename,"afsdn",sector);
%n_fsd = fsd_converter(filename,"afsdn","fsdrad",sector);

raw_data = case1;

alpha = 0.66; % Rothrock and Thorndike (1984)
% Get tracer array
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
afsdn = data_format_sector(filename,"afsdn",sector);
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                trcrn(i,j,k,n) = raw_data(i,j,k,n);%.*floe_binwidth(k)./aicen_data(i,j,n); % something x Length: (something m)
            end
        end
    end
end
% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*floe_rad_c.^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                nfstd(i,j,k,n) = trcrn(i,j,k,n)/floe_area_c(k); % Dimensionless: something / Length
            end
        end
    end
end
% Intergrate nfstd to get nfsd

for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            processed_data(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
        end
    end
end

icemask = aice2 > 0.01;
data2(:,:,:) = processed_data;

for i = 1:Nf
    temp = data1(:,:,i);
    data2_overspace(i) = mean(temp(icemask),'all');
end
% No forcing
filename = 'cases/forcingoff/history/iceh.2005-09-30.nc';
case1 = data_format_sector(filename,"afsdn",sector);
%n_fsd = fsd_converter(filename,"afsdn","fsdrad",sector);

raw_data = case1;

alpha = 0.66; % Rothrock and Thorndike (1984)
% Get tracer array
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
afsdn = data_format_sector(filename,"afsdn",sector);
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                trcrn(i,j,k,n) = raw_data(i,j,k,n);%.*floe_binwidth(k)./aicen_data(i,j,n); % something x Length: (something m)
            end
        end
    end
end
% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*floe_rad_c.^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            for n = 1:Nc
                nfstd(i,j,k,n) = trcrn(i,j,k,n)/floe_area_c(k); % Dimensionless: something / Length
            end
        end
    end
end
% Intergrate nfstd to get nfsd

for i = 1:nx
    for j = 1:ny
        for k = 1:Nf
            processed_data(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
        end
    end
end

icemask = aice2 > 0.01;
data3(:,:,:) = processed_data;

for i = 1:Nf
    temp = data1(:,:,i);
    data3_overspace(i) = mean(temp(icemask),'all');
end% Data reading done

% Plotting


for i = 1
    cat = i;
    nexttile
    %[len,wid] = size(data);
    hold on
    % FIX THIS data_fsd = afsd.ave(:,i).*floe_binwidth'./(4*alpha*NFSD);
    %data_number_fsd = [data.afsd.ave(:,7,1).*floe_binwidth'./(4*alpha*NFSD.^2)];
    p(i) = plot(log10(NFSD),log(data1_overspace*10^6),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
    hold on 
    plot(log10(NFSD),log(data2_overspace*10^6),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
    plot(log10(NFSD),log(data3_overspace*10^6),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
    %data_number_fsd = [data.afsd.ave(:,7,2).*floe_binwidth'./(4*alpha*NFSD.^2)];
    %plot(log10(NFSD),log(data_number_fsd),'-.*','MarkerFaceColor', '#77AC30','LineWidth',1.5);
    plot(log10(NFSD(1:12)),log10(roach_data_sep),'-o','MarkerFaceColor', 'k','Color', 'k','LineWidth',3)
    if i == 1
        legend({'WIM on (Sep)','WIM off (Sep)', "Forcing off (Sep)",'Roach et al. (2018) Sep results'})
    end
    grid on
    %title(sprintf("July",i))
    xticks(log10(xtick))
    xticklabels(xticklab)
    xlabel('Floe radius (m)')
    yticks(log10(ytick))
    yticklabels(yticklab)
    ylim([-30,5])
    %xlim([0,3*10^3]);
end


%xlabel(t5,'FSD cats','FontSize',14)
%ylabel(t5,'Number FSD','FontSize',14)
%title5 = title(t5,strcat(case_name, " Southern hemisphere number FSD, F^N(r)dh"));
%title5.Interpreter = 'Latex';
f5.Position = [1200 100 600 2000];
hold off

%% Figure 2 (h)
f6 = figure(6);
t6 = tiledlayout(2,5);
t6.TileSpacing = 'compact';
alpha = 0.66;
cust_bounds =  max(NFSD);
xtick = 10.^(0:4);
xticklab = cellstr(num2str(round(log10(xtick(:))), '10^{%d}'));
ytick = 10.^(-9:5);
yticklab = cellstr(num2str(round(log10(ytick(:))), '10^{%d}'));

for i = 1:10
    cat = i;
    nexttile
    %[len,wid] = size(data);
    hold on
    % FIX THIS data_fsd = afsd.ave(:,i).*floe_binwidth'./(4*alpha*NFSD);
    p(i) = plot(log10(NFSD),log(data_fsd),'-o','MarkerFaceColor', 'k','Color', 'k');
    grid on
    title(sprintf("Day %d",i))
    xticks(log10(xtick))
    xticklabels(xticklab)
    xlabel('Floe radius (m)')
    yticks(log10(ytick))
    yticklabels(yticklab)
    %xlim([0,3*10^3]);
end

%xlabel(t5,'FSD cats','FontSize',14)
%ylabel(t5,'Number FSD','FontSize',14)
%title5 = title(t5,strcat(case_name, " Southern hemisphere number FSD, F^N(r)dh"));
%title5.Interpreter = 'Latex';
f5.Position = [1200 100 1200 500];
hold off



%% Change in ra


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