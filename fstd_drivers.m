%% What drives the changes in the FSTD?
%% Read in the data.
clear
close all
addpath functions
addpath packages/quiverwcolorbar
clc
% Parameters
sector = "EA";
grid = 'gx1';
case_name = 'ocnatmo';
filedir = '/Volumes/NoahDay5TB/cases/ocnatmo/history/iceh.';
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

initial_date.day = 1;
initial_date.month = 5;
initial_date.year = 2006; % 2005 is a spin-up year
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

[aice, sector_mask] = data_format_sector(filename,"aice",sector);
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

all_cat = 1;

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
