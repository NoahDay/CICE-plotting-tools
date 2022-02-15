%% FSTD_VALIDATION
% This script will validate the FSTD construction through CICE history
% fields to make sure everything is behaving as expected.
% Aim: To get from afsdn to aice
% Method:
% 1. Integrate aicen to aice
% 2. Integreate 


clear all
close all
addpath functions
addpath packages/bedmap2_toolbox_v4
filename = 'cases/fixedwaves/history/iceh.2005-01-01.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';
ssd = 0;
sector = "EA";

%% Preamble
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'fixedwaves';
if ssd == 1
    ssd_dir = '/Volumes/Noah_SSD/run_data';
    filedir = strcat(ssd_dir,case_name);
else
    filedir = strcat('cases/',case_name);
end
grid = 'gx1'; 
time_period = 'd'; %'1','d','m','y'
season = "winter";

day = 1;
month_init = 1;
year = 2005;
date = sprintf('%d-0%d-0%d', year, month_init, day);
months = [1];
dim = 2;
[lat,lon,row] = grid_read(grid);

figcount = 0;

date = sprintf('%d-0%d-0%d', year, month_init, day);
filename = strcat(filedir,"/history/iceh.",date,".nc");
[len, wid] = size(lat);
dim = 3;
variable = "afsd";
NFSD = ncread(filename,"NFSD");
fsd_data = data_format_sector(filename,variable,sector,dim);
%     idx = fsd_data > 10; 
%     fsd_miz = fsd_data;
%     fsd_miz(idx) = 0.0;
for i = 1:len
    for j = 1:wid
        fsd_pdf = zeros(1,length(NFSD));
        for k = 1:length(NFSD)
            fsd_pdf(k) = fsd_data(i,j,k);
        end
        fsd_miz(i,j) = fsd_pdf*(NFSD.^2*pi); % Fraction x Area
        %if ~isnan(fsd_pdf*(NFSD.^2*pi))
        %    temp = fsd_pdf.*(NFSD.^2*pi);
        %end
    end
end
figcount = figcount + 1;
figure(figcount)
map_plot(fsd_miz,variable,sector,grid,[0,1000]);

% Normalize FSD per cell
variable = "afsd";
NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
dim = 3;
fstd_data = data_format_sector(filename,variable,sector,dim); % (x,y,fsd,itd)
for i = 1:len
    for j = 1:wid
        fsd_pdf = zeros(1,length(NFSD));
        for k = 1:length(NFSD)
%                 % Normalize ITD
%                 for l = 1:length(NCAT)
%                     itd_pdf(l) = fstd_data(i,j,k,l)/sum(fstd_data(i,j,k,:));
%                 end
%                 % Integrate the ITD
%                 Int_itd = itd_pdf*NCAT;
%                 fsd_data = fstd_data(:,:,:,1);
            % Intergrate w.r.t. FSD
            fsd_pdf(k) = fsd_data(i,j,k)/sum(fsd_data(i,j,:));
        end
        fsd_miz(i,j) = fsd_pdf*NFSD; % Weighted Average
    end
end
figcount = figcount + 1;
figure(figcount)
map_plot(fsd_miz,variable,sector,grid,[0,1000]);


% FSDrad 
variable = "fsdrad";
fsdrad_data = data_format_sector(filename,variable,sector,2);
figcount = figcount + 1;
figure(figcount)
map_plot(fsdrad_data,variable,sector,grid,[0,1000]);

% AICE 
variable = "aice";
aice_data = data_format_sector(filename,variable,sector,2);
figcount = figcount + 1;
figure(figcount)
map_plot(aice_data,variable,sector,grid,[0,1]);

figcount = figcount + 1;
figure(figcount)
map_plot(fsdrad_data-fsd_miz,variable,sector,grid,[0,1000]);
