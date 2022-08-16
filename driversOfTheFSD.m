% Setup
clear all
clc
close all
addpath functions
addpath /Users/noahday/GitHub/CICE-analyser/processing

historydir = '/Users/noahday/maths1/spectrum_fixed/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end

[lat,lon] = grid_read("gx1");
NFSD = ncread(filenames(1,:),"NFSD");
NCAT = ncread(filenames(1,:),"NCAT");
[floe_binwidth, floe_rad_l, floe_rad_h, floe_area_binwidth] = cice_parameters(NFSD);
%% Read in data
sector = "SA";
close all
clear data
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
for i = 1:n_files
    [dafsd.latm, sector_mask] = data_format_sector(filenames(i,:),"dafsd_latm",sector);
    dafsd.latg = data_format_sector(filenames(i,:),"dafsd_latg",sector);
    dafsd.weld = data_format_sector(filenames(i,:),"dafsd_weld",sector);
    dafsd.newi = data_format_sector(filenames(i,:),"dafsd_newi",sector);
    dafsd.wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
    aice  = data_format_sector(filenames(i,:),"aice",sector);
    afsd = data_format_sector(filenames(i,:),"afsd",sector);

    % Integrate over space
    int.latm(i,:) = integrate_dafsd_hemi(dafsd.latm, aice, floe_binwidth);
    int.latg(i,:) = integrate_dafsd_hemi(dafsd.latg, aice, floe_binwidth);
    int.newi(i,:) = integrate_dafsd_hemi(dafsd.newi, aice, floe_binwidth);
    int.weld(i,:) = integrate_dafsd_hemi(dafsd.weld, aice, floe_binwidth);
    int.wave(i,:) = integrate_dafsd_hemi(dafsd.wave, aice, floe_binwidth);
    data(:,:,i) = [int.latm(i,:); int.latg(i,:); int.newi(i,:); int.weld(i,:); int.wave(i,:)];
    int.afsd(i,:) = integrate_afsd_hemi(afsd, aice, floe_binwidth);
end
%% Plotting change in afsd due to each process [dimless]
close all
colour.latm = [0.8500, 0.3250, 0.0980];
colour.latg = [0.4660, 0.6740, 0.1880];
colour.newi = [0.9290, 0.6940, 0.1250];
colour.weld = [0.3010, 0.7450, 0.9330];
colour.wave = [0, 0.4470, 0.7410];
ymax = max(max(abs(sum(data))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
fig = figure;
conFigure(11)
clear temp
year_idx = 12;
for i = 1:11
    temp(:,:) = data(:,:,i+year_idx);
    [len,wid] = size(temp);
    subplot(3,6,i)
    b = bar(1:length(NFSD),temp,'stacked');
    b(1).FaceColor = colour.latm;
    b(2).FaceColor = colour.latg;
    b(3).FaceColor = colour.newi;
    b(4).FaceColor = colour.weld;
    b(5).FaceColor = colour.wave;
    title(month_label(i))

    ylim([-ymax,ymax])
    
end
i = i + 1;
subplot(3,6,i)
 b = bar(1:length(NFSD),nan(len,wid),'stacked');
    b(1).FaceColor = colour.latm;
    b(2).FaceColor = colour.latg;
    b(3).FaceColor = colour.newi;
    b(4).FaceColor = colour.weld;
    b(5).FaceColor = colour.wave;
axis off
legend('Lat. melt','Lat. growth','New floes','Welding','Breakup','Orientation','vertical','Location','west')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Change in SIC [\%/day]');
xlabel(han,'Floe size categories');
%title(han,'yourTitle');

%% Plotting change in afsd due to each process [1/m]
close all
colour.latm = [0.8500, 0.3250, 0.0980];
colour.latg = [0.4660, 0.6740, 0.1880];
colour.newi = [0.9290, 0.6940, 0.1250];
colour.weld = [0.3010, 0.7450, 0.9330];
colour.wave = [0, 0.4470, 0.7410];
ymax = 0.001;%max(max(abs(sum(data))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
fig = figure;
conFigure(11)
clear temp
year_idx = 24;
for i = 1:12
    temp(:,:) = data(:,:,i+year_idx);
    temp = temp./floe_binwidth;
    % Log the positive values
    %idx = temp > 0;
    %temp(idx) = log10(temp(idx));
    % Log the negative values
    %idx = temp < 0;
    %temp(idx) = -log10(-temp(idx));

    [len,wid] = size(temp);
    subplot(2,7,i)
    b = bar(1:length(NFSD),temp,'stacked');
    %set(gca,'YScale','log')
    b(1).FaceColor = colour.latm;
    b(2).FaceColor = colour.latg;
    b(3).FaceColor = colour.newi;
    b(4).FaceColor = colour.weld;
    b(5).FaceColor = colour.wave;
    title(month_label(i))
    %ytick([-5, -4, -3, -2, -1, 0])
    %yticklabels(['10^(-5)', '10^(-4)', '10^(-3)', '10^(-2)', '10^(-1)', '1'])
    ylim([-ymax,ymax])
    
end
i = i + 1;
subplot(2,7,[i,i+1])
 b = bar(1:length(NFSD),nan(len,wid),'stacked');
    b(1).FaceColor = colour.latm;
    b(2).FaceColor = colour.latg;
    b(3).FaceColor = colour.newi;
    b(4).FaceColor = colour.weld;
    b(5).FaceColor = colour.wave;
axis off
legend('Lat. melt','Lat. growth','New floes','Welding','Breakup','Orientation','vertical','Location','west')

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Change in FSD [m$^{-1}$day$^{-1}$]');
xlabel(han,'Floe size categories');
%title(han,'yourTitle');

%% Plotting change of FSD over time
close all
ymax = max(max(abs(int.afsd(:,:))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
conFigure(14.5,5)
fig = figure;

clear temp
for i = 1:n_files
    temp(:,:) = int.afsd(i,:);
    [len,wid] = size(temp);
    subplot(2,6,i)
    b = bar(1:length(NFSD),temp);
    title(month_label(i))
    xticklabels(round(NFSD))
    ylim([0,ymax])
    
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Fraction of SIC');
xlabel(han,'Floe radius [m]');
%title(han,'yourTitle');
%% Changes in newi over a year
close all
ymax = max(max(abs(int.newi(:,:))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
conFigure(14.5,5)
fig = figure;

clear temp
for i = 1:12 % Month index
    temp(:,:) = int.newi(i,:); % Each month
    [len,wid] = size(temp);
    subplot(2,6,i)
    b = bar(1:length(NFSD),temp);
    b.FaceColor = colour.newi;
    title(month_label(i))
    xticklabels(round(NFSD))
    ylim([-ymax,ymax])
    
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Fraction of SIC');
xlabel(han,'Floe radius [m]');
%title(han,'yourTitle');



%% Changes in welding over a year
close all
ymax = max(max(abs(int.weld(:,:))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
conFigure(14.5,5)
fig = figure;

clear temp
for i = 1:12 % Month index
    temp(:,:) = int.weld(i,:); % Each month
    [len,wid] = size(temp);
    subplot(2,6,i)
    b = bar(1:length(NFSD),temp);
    b.FaceColor = colour.weld;
    title(month_label(i))
    xticklabels(round(NFSD))
    ylim([-ymax,ymax])
    
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Fraction of SIC');
xlabel(han,'Floe radius [m]');
%title(han,'yourTitle');

%% Changes in wave over a year
close all
ymax = max(max(abs(int.wave(:,:))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
conFigure(14.5,5)
fig = figure;

clear temp
for i = 1:12 % Month index
    temp(:,:) = int.wave(i,:); % Each month
    [len,wid] = size(temp);
    subplot(2,6,i)
    b = bar(1:length(NFSD),temp);
    b.FaceColor = colour.wave;
    title(month_label(i))
    xticklabels(round(NFSD))
    ylim([-ymax,ymax])
    
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Fraction of SIC');
xlabel(han,'Floe radius [m]');
%title(han,'yourTitle');



%% Plotting change of first FSD category over the year
close all
ymax = max(max(abs(int.afsd(:,:))))*1.2;
month_label = {'J','F','M','A','M','J','J','A','S','O','N','D','J','F','M','A','M','J','J','A','S','O','N','D'};
colour.latm = [0.8500, 0.3250, 0.0980];
colour.latg = [0.4660, 0.6740, 0.1880];
colour.newi = [0.9290, 0.6940, 0.1250];
colour.weld = [0.3010, 0.7450, 0.9330];
colour.wave = [0, 0.4470, 0.7410];
conFigure(20,5)
fig = figure;
ymax = max(abs(sum(sum(data(:,1,:)))))*1.2;
clear temp
for i = 1:n_files
    temp(:,i) = data(:,1,i);
end
    
b = bar(1:12,temp(:,13:24),"stacked");
b(1).FaceColor = colour.latm;
b(2).FaceColor = colour.latg;
b(3).FaceColor = colour.newi;
b(4).FaceColor = colour.weld;
b(5).FaceColor = colour.wave;
ylim([-ymax,ymax])
xticklabels(month_label)
%legend('Lat. melt','Lat. growth','New floes','Welding','Breakup','Orientation','vertical','Location','southeast')
%set(gca,'XTick',[])

legend('Lat. melt','Lat. growth','New floes','Welding','Breakup','Orientation','vertical','Location','northeast')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Change in SIC [\%/day]');
title('Change in first FSD category')


%%

close all
i = 1
swh = data_format_sector(filenames(i,:),"wave_sig_ht",sector);
figure
[w, a, ~] = map_plot(swh,"wave_sig_ht","SH"); 


swh_mask = swh> 10^(-10);
swh(~swh_mask) = NaN;
figure
[w, a, ~] = map_plot(double(swh_mask),"wave_sig_ht","SH"); 
caxis([0,1])


fsdrad = data_format_sector(filenames(i,:),"fsdrad",sector);
figure
[w, a, ~] = map_plot(fsdrad,"fsdrad","SH"); 
caxis([0,850])

%%
close all
i = 21;
sector = "SH";
conFigure(13)
dafsd_newi = data_format_sector(filenames(i,:),"dafsd_newi",sector);
dafsd_newi12(:,:) = dafsd_newi(:,:,12).*floe_binwidth(12);
figure
[w, a, ~] = map_plot(dafsd_newi12,"dafsd","SH"); 
caxis([-0.05,0.05])
cmocean('balance')

dafsd_wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
dafsd_wave12(:,:) = dafsd_wave(:,:,12).*floe_binwidth(12);
figure
[w, a, ~] = map_plot(dafsd_wave12,"dafsd","SH"); 
caxis([-0.05,0.05])
cmocean('balance')

dafsd_wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
dafsd_wave1(:,:) = dafsd_wave(:,:,1).*floe_binwidth(1);
figure
[w, a, ~] = map_plot(dafsd_wave1,"dafsd","SH"); 
caxis([-0.005,0.005])
cmocean('balance')



fsdrad = data_format_sector(filenames(i,:),"fsdrad","SH");
aice = data_format_sector(filenames(i,:),"aice","SH");
idx = aice > 0.1;
fsdrad(~idx) = NaN;
figure
[w, a, ~] = map_plot(fsdrad.*aice,"fsdrad","SH"); 
caxis([0,850])
%cmocean('balance')



afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd1(:,:) = afsdn(:,:,1,1);
idx = aice > 0.1;
afsd1(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd1.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of smallest floe size and thickness';
exportgraphics(f,'pancakeprop.pdf','ContentType','vector')


afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd2(:,:) = afsdn(:,:,2,1);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(2)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSTD(2,1)';
exportgraphics(f,'fstd2prop.pdf','ContentType','vector')


afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd2(:,:) = afsdn(:,:,1,2);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSTD(1,2)';
exportgraphics(f,'fstd1-2prop.pdf','ContentType','vector')

afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd2(:,:) = afsdn(:,:,1,3);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSTD(1,3)';
exportgraphics(f,'fstd1-3prop.pdf','ContentType','vector')

afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd2(:,:) = afsdn(:,:,1,4);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSTD(1,4)';
exportgraphics(f,'fstd1-4prop.pdf','ContentType','vector')


afsdn = data_format_sector(filenames(i,:),"afsdn","SH");
afsd2(:,:) = afsdn(:,:,1,5);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSTD(1,5)';
exportgraphics(f,'fstd1-5prop.pdf','ContentType','vector')

afsd = data_format_sector(filenames(i,:),"afsd","SH");
afsd1(:,:) = afsd(:,:,1);
idx = aice > 0.1;
afsd1(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd1.*floe_binwidth(1)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of smallest floe size';
exportgraphics(f,'smallprop.pdf','ContentType','vector')


afsd = data_format_sector(filenames(i,:),"afsd","SH");
afsd2(:,:) = afsd(:,:,2);
idx = aice > 0.1;
afsd2(~idx) = NaN;
f = figure;
[w, a, ~] = map_plot(afsd2.*floe_binwidth(2)./aice,"aice","SHplain"); 
caxis([0,1])
a.Label.String = 'Proportion of FSD(2)';
exportgraphics(f,'smallfsd2.pdf','ContentType','vector')
%%
pan_prop = afsd1.*floe_binwidth(1)./aice;

f = figure;
[w, a, ~] = map_plot(afsd1.*floe_binwidth(1)./aice,"aice","SHplain"); 
a.Label.String = 'Proportion of pancake ice';



%% Transect FSD
close
ind_lon = 28;
clear temp afsd
for i = 7
    [dafsd.latm, sector_mask] = data_format_sector(filenames(i,:),"dafsd_latm",sector);
    dafsd.latg = data_format_sector(filenames(i,:),"dafsd_latg",sector);
    dafsd.weld = data_format_sector(filenames(i,:),"dafsd_weld",sector);
    dafsd.newi = data_format_sector(filenames(i,:),"dafsd_newi",sector);
    dafsd.wave = data_format_sector(filenames(i,:),"dafsd_wave",sector);
    aice  = data_format_sector(filenames(i,:),"aice",sector);
    afsd_data = data_format_sector(filenames(i,:),"afsd",sector);

    for nf = 1:length(floe_binwidth)
        temp(:,:) = dafsd.latg(:,:,nf);
        afsd.latg(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.latm(:,:,nf);
        afsd.latm(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.newi(:,:,nf);
        afsd.newi(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.weld(:,:,nf);
        afsd.weld(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
        temp(:,:) = dafsd.wave(:,:,nf);
        afsd.wave(:,:,nf,i) = temp(:,:).*floe_binwidth(nf);
        clear temp
    end
end


[len,wid] = size(aice);
for k = 1
    for i = 1:len
        for j = 1:wid
            afsd.latg2d(i,j,k) = sum(afsd.latg(i,j,:,k));
            afsd.latm2d(i,j,k) = sum(afsd.latm(i,j,:,k));
            afsd.newi2d(i,j,k) = sum(afsd.newi(i,j,:,k));
            afsd.weld2d(i,j,k) = sum(afsd.weld(i,j,:,k));
            afsd.wave2d(i,j,k) = sum(afsd.wave(i,j,:,k));
        end
    end
end
%data = [dafsd.latm(ind_lon,:); dafsd.latg(ind_lon,:); dafsd.newi(ind_lon,:); dafsd.weld(ind_lon,:); dafsd.wave(ind_lon,:)];
data = [afsd_data(ind_lon,:,1).*floe_binwidth(1); afsd_data(ind_lon,:,2).*floe_binwidth(2); afsd_data(ind_lon,:,3).*floe_binwidth(3); afsd_data(ind_lon,:,4).*floe_binwidth(4); afsd_data(ind_lon,:,5).*floe_binwidth(5); afsd_data(ind_lon,:,6).*floe_binwidth(6); afsd_data(ind_lon,:,7).*floe_binwidth(7); afsd_data(ind_lon,:,8).*floe_binwidth(8); afsd_data(ind_lon,:,9).*floe_binwidth(9); afsd_data(ind_lon,:,10).*floe_binwidth(10); afsd_data(ind_lon,:,11).*floe_binwidth(11); afsd_data(ind_lon,:,12).*floe_binwidth(12)];
x = -lat(ind_lon,:);

conFigure(20,2)
figure
area(repmat(x,12,1)',data')
xlim([60,70])
title(sprintf('Transect at %g E', lon(ind_lon,1)))
xlabel('Latitude [$^\circ$ S]')
ylabel('SIC [$\%$]')
legend(["1","2",'3','4','5','6','7','8','9','10','11','12'],'Location','Northwest')


data = [afsd.latg2d(ind_lon,:); afsd.latm2d(ind_lon,:); afsd.newi2d(ind_lon,:); afsd.weld2d(ind_lon,:); afsd.wave2d(ind_lon,:)]; 
conFigure(20,2)
figure
bar(repmat(x,5,1)',data')
xlim([60,70])
title(sprintf('Transect at %g E', lon(ind_lon,1)))
xlabel('Latitude [$^\circ$ S]')
ylabel('SIC [$\%$]')
legend(["Lat. Melt", "Lat. Growth", "New ice", "Weld", "Wave"],'Location','Northwest')

%% Functions
function [data_out] = integrate_dafsd_hemi(dafsd,aice,floe_binwidth)

  % Only take data from that sector
    idx = aice > eps; 
    for nf = 1:length(floe_binwidth)
        temp(:,:) = dafsd(:,:,nf);
        dafsd_tally(:,nf) = temp(idx).*floe_binwidth(nf);
        clear temp
    end
    
    data_out = mean(dafsd_tally);

end

function [data_out] = integrate_afsd_hemi(dafsd,aice,floe_binwidth)

  % Only take data from that sector
    idx = aice > 0.1; 
    for nf = 1:length(floe_binwidth)
        temp(:,:) = dafsd(:,:,nf);
        dafsd_tally(:,nf) = temp(idx).*floe_binwidth(nf)./aice(idx);
        clear temp
    end
    
    data_out = mean(dafsd_tally);

end






function [data_out] = dfsdrad(lat,NFSD,floe_binwidth,afsd)

    [len,wid] = size(lat);
    % Calculate aice
    for i = 1:len
        for j = 1:wid
            temp(:) = afsd(i,j,:);
            fsdrad(i,j) = sum(temp.*floe_binwidth'.*NFSD');
            clear temp
        end
    end
    clear len wid
    
    % Define ice masks
    SIC = eps; % 1% ice mask
    ice_mask = fsdrad > SIC;
    SH_mask = lat < 0;
    
    % Apply ice masks
    temp(:,:) = fsdrad(:,:).*SH_mask;
    %temp(~ice_mask) = NaN;


    data_out = temp;
end



