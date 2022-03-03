% 1. Plot the FSD
%   a) is this a PDF?

% 2. What FSTD categories consitute pancake ice?

%% 1. Preamble
close all
clear all
clc
addpath functions

user = 'noahday'; %a1724548, noahday, Noah
case_name = 'ocntest';%'ocnforcing';
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

NCAT = ncread(filename,"NCAT");
NFSD = ncread(filename,"NFSD");
lon_pos = 28;
lat_pos = edge(lon_pos)-10;
floe_binwidth = [5.2438,8.9763,14.7711,23.3545,35.4569,51.6493,72.1173,96.4015,123.1658,150.0742,173.8638,190.6718,397.7316,479.1093,649.9598,881.7363];
thick_binwidth = NCAT - [0;NCAT(1:end-1)];

%% ITD
% Sum a_{in} = 1
dim = 3;
aicen_data = data_format_sector(filename,"aicen",sector,dim);

Nc = numel(NCAT);

aicen(1:Nc) = aicen_data(lon_pos,lat_pos,:);
conc = aice_data(lon_pos,lat_pos)
a_in = aicen/conc
fprintf('The sum _{n=0}^{N_c} a_{in} is: %g \n', sum(a_in))

%% FSD
dim = 3;
afsd_data = data_format_sector(filename,"afsd",sector,dim);

afsd(1:numel(NFSD)) = afsd_data(lon_pos,lat_pos,:);
% 
fsd = afsd.*floe_binwidth/conc;
%fprintf('The normalised FSD is: %g \n', fsd')

bar(1:16,fsd,'hist')
xlabel('FSD radius')
ylabel('Fraction of sea ice belonging to each category')
Nf = numel(NFSD);
for i = 1:Nf
    lbls(i,:) = num2str(floor(NFSD(i)));
end

xticklabels({lbls})
%% AFSDN
% Sum_{Nc} Sum_{Nf} a_{in} F_{in,k} = 1
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


%% JUNK
sum(afsdn')
afsd
%afsdn_test = afsdn./sum(sum(afsdn))

Nf = numel(NFSD);
for nf = 1:Nf
   test(nf) = sum(afsdn(nf,:)./aicen);
end
test

for nf = 1:numel(NFSD)
    for nc = 1:numel(NCAT)
        trcrn(nf,nc) = afsdn(nf,nc)./(aicen(nc)/floe_binwidth(nf));
    end
end

sum(sum(afsdn))/conc
% for nf = 1:numel(NFSD)
%    ain_calc(nf) = sum(afsdn(:,nf)./(4*0.66*NFSD.^2)); 
% end

% Show that the sum over NCAT and NFSD = 1
% Sum_{nc} Sum_{nf} a_{in} F_{in,k} = 1
% F_{in,k} is the fraction of ice belonging to to thickness category n with
% lateral floe belonging to floe size class k

%sum(afsdn)
%sum(sum(afsdn))
%temp = [NFSD';NFSD';NFSD';NFSD';NFSD'];
%temp2 = [NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT,NCAT, NCAT];
%pdf = afsdn.*temp.*temp2;
%pdf = fsd/aice_data(lon_pos, lat_pos);

%plot(NFSD,pdf)
%sum(sum(pdf))

% surf(afsdn)
% 
% xlabel('ITD')
% ylabel('FSD')
% zlabel('Probability Density')

for nf = 1:Nf
    int_itd = 0;
    for nc = 1:Nc
      int_itd = int_itd +  afsdn(nf,nc);
    end
    pdf(nf) = int_itd/floe_binwidth(nf);
end
pdf
tarea_data = data_format_sector(filename,"tarea",sector,dim);