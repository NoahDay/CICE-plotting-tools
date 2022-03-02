% 1. Plot the FSD
%   a) is this a PDF?

% 2. What FSTD categories consitute pancake ice?

%% 1. Preamble
close all
clear all
addpath functions

user = 'noahday'; %a1724548, noahday, Noah
case_name = 'ocnforcing';
grid = 'gx1'; 
day = 3;
month = 7;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
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

lon_pos = 30;

% Read in data

dim = 3;
afsd_data = data_format_sector(filename,"afsd",sector,dim);
NFSD = ncread(filename,"NFSD");

% Take the FSD at the ice edge (15%)
fsd(1:numel(NFSD)) = afsd_data(lon_pos, edge(lon_pos)-10,:);

sum(fsd)

% AFSDN
dim = 4;
afsdn_data = data_format_sector(filename,"afsdn",sector,dim);
NCAT = ncread(filename,"NCAT");

% Take the FSD at the ice edge (15%)
for nf = 1:numel(NFSD)
    for nc = 1:numel(NCAT)
        afsdn(nc,1:16) = afsdn_data(lon_pos, edge(lon_pos)-10,:,nc);
    end
end

sum(afsdn)
sum(sum(afsdn))
temp = [NFSD';NFSD';NFSD';NFSD';NFSD'];
temp2 = [NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT];
pdf = afsdn.*temp.*temp2;
plot(NFSD,pdf)
sum(sum(pdf))