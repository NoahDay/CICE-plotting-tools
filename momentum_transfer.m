%% Calculate the momentum of the ice from the components and equation

%% Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar

% Parameters
sector = "SA";
grid = 'gx1';
case_name = 'momentum';
filedir = 'cases/momentum/history/iceh.';
%'cases/momentum/history/iceh.'; %'/Volumes/NoahDay5TB/cases/momentum/history/iceh.2009-09-30.nc';
[lat,lon,row] = grid_read(grid);

% Make sector mask
%[len,wid] = size(lat);
fig_count = 0;


coords = sector_coords(sector);
% Take the location to be the centre of the sector
sector = "SA";%struct('coords',sector_coords(sector),'centre_lon',(coords(1,1)+coords(2,1))/2,'north_lat',[],'south_lat',[]);
clear coords

initial_date.day = 3;
initial_date.month = 7;
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

lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);

SIC = 0.15;
[lat_ice_edge, lon_ice_edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);

% 1. Ice drift
% Variables: uvel, vvel

uvel_data = data_format_sector(filename,"uvel",sector);
vvel_data = data_format_sector(filename,"vvel",sector);

uvel_vec = reshape(uvel_data,1,[]);
vvel_vec = reshape(vvel_data,1,[]);
ice_vel_mag = sqrt(uvel_data.^2 + vvel_data.^2);



t = tiledlayout(1,5);



%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.ice = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,uvel_vec,vvel_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Ice drift velocity - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)
% 2. Wind stress on ice
% Variables: strairx, strairy

% Plot a map of the stress vectors over the top of aice

strairx_data = data_format_sector(filename,"strairx",sector);
strairy_data = data_format_sector(filename,"strairy",sector);

strairx_vec = reshape(strairx_data,1,[]);
strairy_vec = reshape(strairy_data,1,[]);
strair_mag = sqrt(strairx_data.^2 + strairy_data.^2);

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.wind = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,strairx_vec,strairy_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Wind stress on ice - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)

% 3. Ocean stress
% Variables: strocnx strocny

strocnx_data = data_format_sector(filename,"strocnx",sector);
strocny_data = data_format_sector(filename,"strocny",sector);

strocnx_vec = reshape(strocnx_data,1,[]);
strocny_vec = reshape(strocny_data,1,[]);
strocn_mag = sqrt(strocnx_data.^2 + strocny_data.^2);

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.ocn = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,strocnx_vec,strocny_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Ocean stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)

% 4a. Internal stress
% Variables: sig1, sig2, sigP

nexttile
strength_data = data_format_sector(filename,"sigP",sector);
figs.cor = map_plot(strength_data,"sigP",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
   % quiverm(lat_vec,lon_vec,strcorx_vec,strcory_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Coriolis stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)


% 4b. Divergence of internal
% Variables: divu

% 5. Coriolis effect
% Variables: strcorx, strcory

strcorx_data = data_format_sector(filename,"strcorx",sector);
strcory_data = data_format_sector(filename,"strcory",sector);

strcorx_vec = reshape(strcorx_data,1,[]);
strcory_vec = reshape(strcory_data,1,[]);
strcor_mag = sqrt(strcorx_data.^2 + strcory_data.^2);

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.cor = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,strcorx_vec,strcory_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Coriolis stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)


% 6. Sea surface slope
% Variables: strtltx, strtlty

t2 = tiledlayout(1,3);
strtltx_data = data_format_sector(filename,"strtltx",sector);
strtlty_data = data_format_sector(filename,"strtlty",sector);

strtltx_vec = reshape(strtltx_data,1,[]);
strtlty_vec = reshape(strtlty_data,1,[]);
strtlt_mag = sqrt(strtltx_data.^2 + strtlty_data.^2);

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.tlt = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,strtltx_vec,strtlty_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Sea surface slope stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)
% 7. Internal ice stress
% Variables: strintx, strinty

strintx_data = data_format_sector(filename,"strintx",sector);
strinty_data = data_format_sector(filename,"strinty",sector);

strintx_vec = reshape(strintx_data,1,[]);
strinty_vec = reshape(strinty_data,1,[]);
strint_mag = sqrt(strintx_data.^2 + strinty_data.^2);

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.int = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,strintx_vec,strinty_vec,pram.quiver_color)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Internal ice stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)


% 8. Sum the stresses
% Variables x: strairx, strocnx, strcorx, strtltx, strintx
% Variables y: strairy, strocny, strcory, strtlty, strinty

sum_vec_x = strairx_vec + strocnx_vec + strcorx_vec + strtltx_vec + strintx_vec;
sum_vec_y = strairy_vec + strocny_vec + strcory_vec + strtlty_vec + strinty_vec;
sum_mag = sqrt(sum_vec_x.^2 + sum_vec_y.^2);
scale = 50;

%fig_count = fig_count + 1;
%figure(fig_count)
nexttile
figs.sum = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    quiverm(lat_vec,lon_vec,sum_vec_x,sum_vec_y)
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("Sum of ice stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    cmocean(pram.colormap,10)

%% All stresses    
fig_count = fig_count + 1;
figure(fig_count)
fig.all = map_plot(aice_data,"aice",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    q1 = quiverm(lat_vec,lon_vec,strairx_vec,strairy_vec,'b');
    q2 = quiverm(lat_vec,lon_vec,strocnx_vec,strocny_vec,'r');
    q3 = quiverm(lat_vec,lon_vec,strcorx_vec,strcory_vec,'g');
    q4 = quiverm(lat_vec,lon_vec,strtltx_vec,strtlty_vec,'m');
    q5 = quiverm(lat_vec,lon_vec,strintx_vec,strinty_vec,'k');
    %q6 = quiverm(lat_vec,lon_vec,uvel_vec,vvel_vec,	'k');
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("All of ice stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    legend([q1(1), q2(1), q3(1), q4(1), q5(1)], 'Air', 'Ocean','Coriolis','Sea surface slope','Internal')
    cmocean(pram.colormap,10)

% FSD radius
fsd_data = data_format_sector(filename,"fsdrad",sector);
idx = fsd_data < eps;
fsd_data(idx) = NaN;
fig_count = fig_count + 1;
figure(fig_count)
fig.all = map_plot(fsd_data,"fsdrad",sector);  
    plotm(lat_ice_edge,lon_ice_edge,'-','color',pram.ice_edge_color,'LineWidth',2)
    q1 = quiverm(lat_vec,lon_vec,strairx_vec,strairy_vec,'b');
    q2 = quiverm(lat_vec,lon_vec,strocnx_vec,strocny_vec,'r');
    q3 = quiverm(lat_vec,lon_vec,strcorx_vec,strcory_vec,'g');
    q4 = quiverm(lat_vec,lon_vec,strtltx_vec,strtlty_vec,'m');
    q5 = quiverm(lat_vec,lon_vec,strintx_vec,strinty_vec,'k');
    %q6 = quiverm(lat_vec,lon_vec,uvel_vec,vvel_vec,	'k');
    t = textm(lat_ice_edge(20),lon_ice_edge(18),sprintf('SIC = %g\n ice edge',pram.min_SIC),'HorizontalAlignment','right');
    t.Color = pram.label_color;
    set(gcf,'Position',[1200 1000 300 400])
    title(strcat("All of ice stress - ",initial_date.char), 'interpreter','latex','FontSize', 14)
    legend([q1(1), q2(1), q3(1), q4(1), q5(1)], 'Air', 'Ocean','Coriolis','Sea surface slope','Internal')
    cmocean(pram.colormap,30)
    %cmap(parula(20))
    caxis([0 1000])
    %set(gca,'ColorScale','log')

%% Stress comparison

initial_date.day = 28;
initial_date.month = 6;
initial_date.year = 2009; 
if initial_date.day > 9
    initial_date.char = sprintf('%d-0%d-%d', initial_date.year, initial_date.month, initial_date.day);
else
    initial_date.char = sprintf('%d-0%d-0%d', initial_date.year, initial_date.month, initial_date.day);
end
datapoints = 10;

filedir = 'cases/momentum';%'/Volumes/NoahDay5TB/cases/momentum';
% 1. Pick a cell along the ice edge

% For now just pick (-55 S, 30 E)
cell_coords = [-68,30];
[cell_lat,cell_lon] = lat_lon_finder(cell_coords(1),cell_coords(2),lat,lon);

% 2. Store all the stress components over time
date = initial_date.char;
for i = 1:datapoints

    filename = strcat(filedir,"/history/iceh.",date,".nc");
    
    % Find ice edge
    aice_data = data_format_sector(filename,"aice",sector);
    [lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);
    %lon_ice_edge(33) = 30
    [cell_lat,cell_lon] = lat_lon_finder(lat_ice_edge(33),lon_ice_edge(33),lat,lon);
    cell_lat = cell_lat-1;

    % Air stresses
    strairx_data = data_format_sector(filename,"strairx",sector);
    strairy_data = data_format_sector(filename,"strairy",sector);
    strair_mag = sqrt(strairx_data.^2 + strairy_data.^2);
    strair_cell(i,:) = [strairx_data(cell_lon,cell_lat), strairy_data(cell_lon,cell_lat)];
    strair_mag_cell(i) = strair_mag(cell_lon,cell_lat);

    % Ocean stresses
    strocnx_data = data_format_sector(filename,"strocnx",sector);
    strocny_data = data_format_sector(filename,"strocny",sector);
    strocn_mag = sqrt(strocnx_data.^2 + strocny_data.^2);
    strocn_cell(i,:) = [strocnx_data(cell_lon,cell_lat), strocny_data(cell_lon,cell_lat)];
    strocn_mag_cell(i) = strocn_mag(cell_lon,cell_lat);

    % Coriolis stresses
    strcorx_data = data_format_sector(filename,"strcorx",sector);
    strcory_data = data_format_sector(filename,"strcory",sector);
    strcor_mag = sqrt(strcorx_data.^2 + strcory_data.^2);
    strcor_cell(i,:) = [strcorx_data(cell_lon,cell_lat), strcory_data(cell_lon,cell_lat)];
    strcor_mag_cell(i) = strcor_mag(cell_lon,cell_lat);
    
    % Sea surface slope stresses
    strtltx_data = data_format_sector(filename,"strtltx",sector);
    strtlty_data = data_format_sector(filename,"strtlty",sector);
    strtlt_mag = sqrt(strtltx_data.^2 + strtlty_data.^2);
    strtlt_cell(i,:) = [strtltx_data(cell_lon,cell_lat), strtlty_data(cell_lon,cell_lat)];
    strtlt_mag_cell(i) = strtlt_mag(cell_lon,cell_lat);
    
    % Internal stresses
    strintx_data = data_format_sector(filename,"strintx",sector);
    strinty_data = data_format_sector(filename,"strinty",sector);
    strint_mag = sqrt(strintx_data.^2 + strinty_data.^2);
    strint_cell(i,:) = [strintx_data(cell_lon,cell_lat), strinty_data(cell_lon,cell_lat)];
    strint_mag_cell(i) = strint_mag(cell_lon,cell_lat);
    
    % Update date
    date = update_date(date);
end

% 3. Plot a time series of the magnitude of stresses over time

end_date = date;

ts.init_date = datetime(initial_date.year,initial_date.month,initial_date.day);
ts.end_date = datetime(str2num(end_date(1:4)),str2num(end_date(6:7)),str2num(end_date(9:10)));
ts.dates = datevec(ts.init_date:ts.end_date);

% All stresses

data_all_stress = [strocn_mag_cell; strcor_mag_cell;strtlt_mag_cell;strint_mag_cell];
ts_str = timeseries_plot(data_all_stress','Stress induced on sea ice','days',char(ts.init_date));
ts_atm = timeseries_plot(strair_mag_cell','Stress induced on sea ice','days',char(ts.init_date));
fig_count = fig_count + 1;
figure(fig_count)
plot(ts_str,'-', 'LineWidth',2)
    set(gcf,'Position',[1200 1000 500 200])
    ylabel('Stress')
    hold on
    plot(ts_atm,'--', 'LineWidth',2)
    legend({'Ocean','Coriolis','Sea surface slope','Internal','Air'},'Location','northwest')
    hold off


% 4. Plot a time series of the direction of stresses over time

data_all_stress_dir = (360/pi)*[atan(strocn_cell(:,2)./strocn_cell(:,1)), atan(strcor_cell(:,2)./strcor_cell(:,1)), atan(strtlt_cell(:,2)./strtlt_cell(:,1)), atan(strint_cell(:,2)./strint_cell(:,1))];
ts_str_dir = timeseries_plot(data_all_stress_dir,'Direction of stress','days',char(ts.init_date));
atm_stress_dir = (360/pi)*atan(strair_cell(:,2)./strair_cell(:,1));
ts_atm_dir = timeseries_plot(atm_stress_dir,'Direction of stress','days',char(ts.init_date));
fig_count = fig_count + 1;
figure(fig_count)
plot(ts_str_dir,'LineWidth',2)
    set(gcf,'Position',[1200 1000 500 200])
    ylabel('Degrees')
    hold on
    plot(ts_atm_dir,'--', 'LineWidth',2)
    legend({'Ocean','Coriolis','Sea surface slope','Internal','Air'},'Location','northwest')
    hold off


%% Functions

function ts = timeseries_plot(data,name,units,startdate)
    ts = timeseries(data,1:length(data));
    ts.Name = name;
    ts.TimeInfo.Units = units;
    ts.TimeInfo.StartDate = startdate;     % Set start date.
    ts.TimeInfo.Format = 'dd mmm yy';       % Set format for display on x-axis.
end

