% 1. Plot the FSD
%   a) is this a PDF?

% 2. What FSTD categories consitute pancake ice?

%% 1. Preamble
close all
clear all
clc
addpath functions

user = 'a1724548'; %a1724548, noahday, Noah
case_name = 'ocntest';%'ocnforcing';
grid = 'gx1'; 
day = 3;
month = 12;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end

dim = 2;
[lat,lon,row] = grid_read(grid);


ssd = 1;
if ssd == 1
    filename = strcat('/Volumes/NoahDay5TB/cases/ocnforcing/history/iceh.',date,".nc");
else
    filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
end

aice_data = data_format_sector(filename,"aice",sector);
lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);
SIC = 0.15;
[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);

NCAT = ncread(filename,"NCAT");
NFSD = ncread(filename,"NFSD");
Nf = numel(NFSD);
lon_pos = 28;
lat_pos = edge(lon_pos)-8;
floe_binwidth = [5.2438,8.9763,14.7711,23.3545,35.4569,51.6493,72.1173,96.4015,123.1658,150.0742,173.8638,190.6718,397.7316,479.1093,649.9598,881.7363];
thick_binwidth = NCAT - [0;NCAT(1:end-1)];

[nx,ny] = size(aice_data);
% ITD
% Sum a_{in} = 1
dim = 3;
aicen_data = data_format_sector(filename,"aicen",sector,dim);

Nc = numel(NCAT);

aicen(1:Nc) = aicen_data(lon_pos,lat_pos,:);
conc = aice_data(lon_pos,lat_pos);
a_in = aicen/conc;
fprintf('The sum _{n=0}^{N_c} a_{in} is: %g \n', sum(a_in))

% FSD

% floe_rad_l = lims(1:nfsd  )
% floe_rad_h = lims(2:nfsd+1)
% floe_binwidth = floe_rad_h - floe_rad_l
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
f_bin_width = floe_rad_h - floe_rad_l;

dim = 3;
afsd_data = data_format_sector(filename,"afsd",sector,dim);

afsd(1:numel(NFSD)) = afsd_data(lon_pos,lat_pos,:);
fsd = afsd.*floe_binwidth/conc;



% AFSDN
% Sum_{Nc} Sum_{Nf} a_{in} F_{in,k} = 1
clc
dim = 4;
afsdn_data = data_format_sector(filename,"afsdn",sector,dim);

% Take the FSD at the ice edge (15%)

for nf = 1:numel(NFSD)
    for nc = 1:numel(NCAT)
        afsdn(nf,nc) = afsdn_data(lon_pos, lat_pos,nf,nc);
    end
end

afsd_cal = sum(afsdn');
afsd_cal.*floe_binwidth; % = fsd

for nc = 1:Nc
    afsdn_norm(:,nc) = afsdn(:,nc).*floe_binwidth';
end


[converted] = fsd_converter(filename,"afsdn","afsd");

fsd_conv = fsd_converter(filename,"afsdn","fsd");
for i = 1:Nf
    fsd_cell(i) = fsd_conv(lon_pos,lat_pos,i);
end
norm_fsd_cell = fsd_cell./aice_data(lon_pos,lat_pos);
scaled = norm_fsd_cell.*(NFSD)';
mode(scaled);
fsdrad = data_format_sector(filename,"fsdrad",sector,2);
fsdrad(lon_pos,lat_pos)

N_fsd = fsd_converter(filename,"afsdn","number_fstd");
N_fsd(lon_pos,lat_pos,:)
% 
% for i = 1:nx
%     for j = 1:ny
%         for k = 1:Nf
%             for n = 1:Nc
%                 afsd_scaled_thickness(i,j,k) = afsdn_data(i,j,k,n)*
%% Plotting

% FSD histogram
figure(1)
bar(1:Nf,fsd,'hist')
xlabel('FSD radius (m)')
ylabel('Fraction of sea ice area')
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
    title(sprintf("FSD at (%g S, %g E)", lat(lon_pos,lat_pos),lon(lon_pos,lat_pos)))
    set(gcf,'Position',[1300 1000 800 400])
    xtickangle(45)
    
        
ncat_vec = [0; NCAT(1:Nc-1)];  
r_NCAT = round(round(NCAT,2),2);
r1_NCAT = round(round(ncat_vec,2),2);
for i = 1:Nc
    %lab2{i} = sprintf('> %g',ncat_vec(i));
    lab2{i} = sprintf('[%g, %g]',r1_NCAT(i),r_NCAT(i));
end
% FSTD
% figure(2)
% surf(afsdn)
%     xlabel('Ice thickness (m)')
%     xticks(1:Nc)
%     xticklabels(lab2)
%     ylabel('FSD radius (m)')
%     yticks(1:Nf)
%     yticklabels(lab)
%     title(sprintf("FSTD at (%g S, %g E)", lat(lon_pos,lat_pos),lon(lon_pos,lat_pos)))

figure(3)
bar3(afsdn_norm)
    xlabel('Ice thickness (m)')
    xticks(1:Nc)
    xticklabels(lab2)
    ylabel('FSD radius (m)')
    yticks(1:Nf)
    yticklabels(lab)
    title(sprintf("FSTD at (%g S, %g E)", lat(lon_pos,lat_pos),lon(lon_pos,lat_pos)))
    set(gcf,'Position',[1300 1000 600 600])
    
figure(4)
    bar(afsdn_norm,'stacked')
    set(gcf,'Position',[1300 1000 800 400])
    [hleg,icons,plots] = legend('show');
    [leg,att] = legendflex(gca, lab2, 'title', 'ITD (m)');
    set(findall(leg, 'string', 'ITD (m)'), 'fontweight', 'bold');
    leg.Title.Visible = 'on';
    xlabel('FSD radius (m)')
    xticks(1:Nf)
    xticklabels(lab)
    xtickangle(45)
    title(sprintf("FSTD at (%g S, %g E) %g cells south of the ice edge", lat(lon_pos,lat_pos),lon(lon_pos,lat_pos), edge(lon_pos) - lat_pos -1))
    ylabel('Fraction of sea ice area')
