% 1D
%clear all
close all
%lon_pos = 30;
clear Hs Hs_nospread Hs_scattering Hs_adj

%lon_required = 90;
%cell_lon_idx = near1(lon(:,1),lon_required); % The cell with longitude closest to lon_required

% Parameters ----------------------------------------------------------------------------------------
    tolice = 1.0e-2;               % conc<tolice treated as zero ice
    tolh   = 1.0e-1;               % h<tolh      treated as zero ice
    toli   = 1.0e-16;              % threshold for the propagation of waves
    gravity = 9.80665;
    puny = 1.0e-11;       


length_transect = 10;

Hs_init = 6;
Tp_init = 10;

dir_spread = 1; % = 1, we are reading in CAWCR data

% 1. Intialise wave spectrum at wave_lat
nw = 31;
clear S_attn S_attn_nospread int_D S_attn_scattering
[S_init,omega,T] = SDF_Bretschneider(Hs_init,Tp_init,nw); 
wavefreq = omega/(2*pi);
dwavefreq = wavefreq(1:end) - [0,wavefreq(1:end-1)];
mwd = 0; % Mean wave direction, rad
thn = 31; % number of theta bins
n = 2.5; % Cosine exponent

if dir_spread == 1
    [S_spread,D,theta_vec] = cosine_spreader(S_init,mwd,thn,n);
    dtheta = theta_vec(1:end) - [0,theta_vec(1:end-1)];
    S_attn(1,:) = S_spread;
    S_attn_nospread(1,:) = S_init;
    S_attn_scattering(1,:) = S_spread;
    S_attn_adj = S_spread;
else
    S_attn(1,:) = S_init;
end

% 2. Propagate the waves, sub_uncoupled          
% increment_floe
%! !DESCRIPTION:
%!
%!  Increase ice floe tracer by scaled timestep length.

 % Length of cell, converting km to m

dwavew = dwavefreq.*(2*pi);
scatter_loc = 10;
n_points = 100;
dist = linspace(0,150,n_points); % In km
delta_dist = [dist - [0 dist(1:end-1)]];
%floe_size_vec = exp(0.05*(dist));

       
        vert_translation = (850+1)/2;
        horz_translation = 50;
        amplitude = (850-1)/2;
        
%for lp_i=1:n_points
%    floe_size_vec(lp_i) = amplitude*tanh(0.001*(lp_i-horz_translation)) + vert_translation;
%end 
%floe_size_vec = 850./(1 + exp(-0.001*x_vec));

floe_size_vec = 850./(1 + exp(-0.1*dist+10));
aice_vec = 0.95./(1 + exp(-0.1*dist+5));
%figure
%plot(dist,floe_size_vec)

%

flg_amp_drop = zeros(1,length(omega));
for i = 1:n_points
    conc = aice_vec(i);
    hice = 1.0;
    floe_size = floe_size_vec(i);% Radius, m
    Lcell = delta_dist(i).*1000; % initialise propagation length

    if conc < tolice 
        % As there is no ice, apply no attenuation
        S_attn(i+1,:) = S_attn(i,:);
    elseif conc > tolice % Attenuate
        % Amplitude drop check
        for om_i = 1:length(omega)
            if flg_amp_drop(om_i) == 0
                T = (2*pi)./omega(om_i);
                lambda = (gravity*(T.^2))/(2*pi); % Wave lengths
                if conc > 0.1 && lambda <= 2*floe_size
                    % Unbroken ice, apply the amplitude drop-off
                    if flg_amp_drop(om_i) == 0
                        disp('Scatter')
                        disp(lambda)
                        disp(floe_size)
                        %[S_attn(i,:) ] = amplitude_dropoff(S_attn(i,:),omega,floe_size);
                        S_attn(i,om_i) = 0.5*S_attn(i,om_i);
                        flg_amp_drop(om_i) = 1;
                        scatter_loc(om_i) = i;
                    end
                    %S_attn(i+1,:) = wave_attenuation(conc,Lcell,nw,S_attn(i,:),omega,"adjusted",floe_size);  
                end
            else
                % Broken ice, no drop-off
                S_attn(i,om_i) = S_attn(i,om_i);
            end
                    
        end    
        S_attn(i+1,:) = wave_attenuation(conc,Lcell,nw,S_attn(i,:),omega,"MBK",floe_size);
    end
    
    Hs(i) = 4*sqrt(sum(S_attn(i,:).*dwavew));

end % i, prop_length


addpath functions
conFigure(11)
figure
yyaxis left
plot(dist,Hs)
ylabel('$H_s$ [m]')
ylim([0,8])
yyaxis right
plot(dist,floe_size_vec./1000)
hold on
plot(dist,aice_vec)
ylabel('$r_a$ [km] \& SIC [-]')
%ylim([0,1000])
xline(dist(scatter_loc),'--',{'$r_a = 100$'},'Interpreter','latex')
set(gca,'YScale','linear')
xlabel('Distance [km]')


 T = (2*pi)./omega;
lambda = (gravity*(T.^2))/(2*pi); % Wave lengths

close all
conFigure(10,1.5)
f = figure;
yyaxis left
plot(dist,Hs)
ylabel('$H_s$ [m]')
hold on
yyaxis right
plot(dist,aice_vec)
plot(dist,floe_size_vec./1000)
ylabel('$r_a$ [km] \& SIC [-]')
for i = 2:9%length(scatter_loc)
    %xline(dist(scatter_loc(i)),'--',{strcat('$\lambda$' ,sprintf('= %g',round(lambda(i),-1)))},'Interpreter','latex','HandleVisibility','off')
end
legend({'$H_s$','SIC','Floe radius ($r_a$)'},'AutoUpdate','off','Location','northoutside','Orientation','horizontal')
set(gca,'YScale','linear')
ylim([0,1])
%ylabel('Normalised values')

xlabel('Distance [km]')
exportgraphics(f,'amp_dropoff.pdf','ContentType','vector')
%addpath /Users/noahday/Documents/MATLAB/matlab2tikz/src/
%matlab2tikz('1d_amp_drop.tex', 'standalone', true);

%%
close all



% Create a 3x4 array of sample data in the range of 0-255.

data = repmat(dist',[1,nw-1]);
data(1:10,:) = 10;
data = zeros(length(dist),nw-1);
for j = 1:nw-1
    for i = 1:length(dist)
        if i >= scatter_loc(j+1)
          
            % Consolidated ice
            data(i,j) = 1;
        end
    end
end


X = repmat(lambda(2:end),[length(dist),1]);
Y = repmat(dist,[1,nw-1]);
conFigure(14,2)
f = figure;
pcolor(data')
%hold on
%plot(1:length(dist),floe_size_vec)
% Initialize a color map array of 256 colors.


yticks(1:5:nw-1)
yticklabels(round(lambda(2:5:end),2,"significant"))
ylabel('$\lambda$ [m]')

xticks(1:10:n_points)
xlabel('Distance from ice edge [km]')
xticklabels(round(dist(1:10:end),0))

colorMap = jet(2);
 cbh = colorbar ; %Create Colorbar
 cbh.Ticks = [0,1] ; %Create 8 ticks from zero to 1
 cbh.TickLabels = {'Unconsolidated','Consolidated'} ;    %Replace the labels of these 8 ticks with the numbers 1 to 8


%%
close all
figure
C = hadamard(20);
pcolor(data)
colormap(gray(3))
axis ij
axis square


%%
close all
floe_size = 10;
[S_test beta_N] = amplitude_dropoff(ones(1,31),omega,floe_size);

figure
plot(omega,S_test)
xlabel('$\omega$')
ylim([-0.2,1.2])
ylabel('$S_{drop}(\omega)/S_{0}(\omega)$')


figure
stairs(omega,beta_N,'linewidth',3)
xlabel('$\omega$')
ylim([-0.1,1.1])
ylabel('$\beta_N$')

gravity = 9.80665;
T = (2*pi)./omega;
lambda = (gravity*(T.^2))/(2*pi); % Wave lengths

idx = 2*floe_size < (lambda/2);
pos = find(idx);

figure
plot(lambda,S_test,'linewidth',3)
set(gca,'XScale','log')
xlabel('$\lambda$')
ylim([-0.1,1.1])
xlim([0,10^3])
xline(lambda(pos(end)),'--',{'$d=\lambda/2$'},'Interpreter','latex','LabelVerticalAlignment','bottom')
ylabel('$S_{drop}(\lambda)/S_{0}(\lambda)$')
title(sprintf('Transmission, $d$ = %g m',floe_size))

%% CICE propagation

% Read data in ----------------------------------------------------------------------------------------
addpath functions
filename_ww3 = '/Users/noahday/GitHub/cice-dirs/input/CICE_data/forcing/access-om2_1deg/CAWCR/MONTHLY/2017/ww3_om2_1deg_201707.nc';
filename = '/Volumes/NoahDay5TB/WIMonAlessandroRun/history/iceh.2017-07-01.nc';
sector = "SH";
grid = 'om2';
aicen_data = data_format_sector(filename,'aicen',sector);
aice_data = sum(aicen_data,3);

swh_data = data_format_cawcr(filename_ww3,'hs');
swh_data = swh_data(:,:,1);
fp_data = data_format_cawcr(filename_ww3,'fp');
fp_data = fp_data(:,:,1);
idx = fp_data < 0;
fp_data(idx) = NaN;
ppd_data = 1./fp_data;
mwd_data = data_format_cawcr(filename_ww3,'dir');
mwd_data = mwd_data(:,:,1);
hte_data = data_format_sector(filename,'HTE',sector);
hice_data = data_format_sector(filename,'hi',sector);
fsdrad_data = data_format_sector(filename,'fsdrad',sector);
[lat,lon] = grid_read(grid);


%lon_pos = 30;
clear Hs Hs_dropoff
Hs = zeros(size(lon));
Hs_dropoff = zeros(size(lon));

%lon_required = 90;
%cell_lon_idx = near1(lon(:,1),lon_required); % The cell with longitude closest to lon_required

% Parameters ----------------------------------------------------------------------------------------
    tolice = 1.0e-2;               % conc<tolice treated as zero ice
    tolh   = 1.0e-1;               % h<tolh      treated as zero ice
    toli   = 1.0e-16;              % threshold for the propagation of waves
    gravity = 9.80665;
    puny = 1.0e-11;       




lon_range = 1:360;
for lon_pos = lon_range
    lat_range = 1:60;
    swh_data_vec= swh_data(lon_pos,lat_range);
    aice_data_vec= aice_data(lon_pos,lat_range);
    hice_data_vec= hice_data(lon_pos,lat_range);
    lat_vec = lat(lon_pos,lat_range);
    

    
    cell_lon_idx = lon_pos;
    aice_gt_puny = find(aice_data_vec>eps);
    if isempty(aice_gt_puny)
        wave_lat = 3;
    else
        wave_lat =  aice_gt_puny(end)+1; % The last cell with aice > puny
    end
    length_transect = wave_lat - 1;
    
    Hs_init = swh_data(cell_lon_idx,wave_lat);
    Tp_init = ppd_data(cell_lon_idx,wave_lat);
    
    dir_spread = 1; % = 1, we are reading in CAWCR data
    
    % 1. Intialise wave spectrum at wave_lat
    nw = 31;
    clear S_attn S_attn_dropoff int_D
    [S_init,omega,T] = SDF_Bretschneider(Hs_init,Tp_init,nw); 
    wavefreq = omega/(2*pi);
    dwavefreq = wavefreq(1:end) - [0,wavefreq(1:end-1)];
    mwd = mwd_data(cell_lon_idx, wave_lat)*(pi/180); % Mean wave direction, rad
    thn = 31; % number of theta bins
    n = 2.5; % Cosine exponent
    
    if dir_spread == 1
        [S_spread,D,theta_vec] = cosine_spreader(S_init,mwd,thn,n);
        dtheta = theta_vec(1:end) - [0,theta_vec(1:end-1)];
        S_attn(1,:) = S_spread;
        S_attn_dropoff(1,:) = S_spread;
    else
        S_attn(1,:) = S_init;
    end
    
    % 2. Propagate the waves, sub_uncoupled          
    % increment_floe
    %! !DESCRIPTION:
    %!
    %!  Increase ice floe tracer by scaled timestep length.
    
    L_data = hte_data(cell_lon_idx,wave_lat:-1:wave_lat-length_transect); % Length of cell, converting km to m
    
    
    dwavew = dwavefreq.*(2*pi);
    conc_vec = aice_data_vec(wave_lat-1:-1:1);
    scatter_pos = find(conc_vec > 0.1);
    if isempty(scatter_pos)
        scatter_loc = 0;
    else
        scatter_loc = scatter_pos(1);
        scatter_loc_vec(lon_pos) = scatter_loc;
    end
    
    prop_length = find(~isnan(aice_data_vec(wave_lat:-1:1)));
    if isempty(prop_length)
        prop_length = 1;
    else
        prop_length = prop_length(end);% Propagation length, cells
        prop_length = wave_lat - 2;
    end

    flg_amp_drop = zeros(1,length(omega));
    for i = 1:prop_length
        conc = aice_data_vec(wave_lat-i);
        hice = hice_data_vec(wave_lat-i);
        floe_size = fsdrad_data(cell_lon_idx,wave_lat-i-1);% Radius, m
        Lcell = L_data(wave_lat-i-1); % initialise propagation length
    
        if conc < tolice 
            % As there is no ice, apply no attenuation
            S_attn(i+1,:) = S_attn(i,:);
            S_attn_dropoff(i+1,:) = S_attn_dropoff(i,:);
        elseif conc > tolice % Attenuate
            % MBK attenuation
            S_attn(i+1,:) = wave_attenuation(conc,Lcell,nw,S_attn(i,:),omega,"MBK",floe_size);
          
            % Amplitude drop check
            for om_i = 1:length(omega)
                if flg_amp_drop(om_i) == 0
                    T = (2*pi)./omega(om_i);
                    lambda = (gravity*(T.^2))/(2*pi); % Wave lengths
                    if conc > 0.1 && lambda <= 2*floe_size
                        % Unbroken ice, apply the amplitude drop-off
                        if flg_amp_drop(om_i) == 0
                            %disp('Scatter')
                            %disp(lambda)
                            %disp(floe_size)
                            %[S_attn(i,:) ] = amplitude_dropoff(S_attn(i,:),omega,floe_size);
                            S_attn_dropoff(i,om_i) = 0.5*S_attn_dropoff(i,om_i);
                            flg_amp_drop(om_i) = 1;
                            scatter_loc(om_i) = i;
                        end  
                    end
                else
                    % Broken ice, no drop-off
                    S_attn_dropoff(i,om_i) = S_attn_dropoff(i,om_i);
                end
                        
            end    
            S_attn_dropoff(i+1,:) = wave_attenuation(conc,Lcell,nw,S_attn_dropoff(i,:),omega,"MBK",floe_size);
        elseif isnan(conc)
            S_attn(i+1,:) = NaN;
            S_attn_dropoff(i+1,:) = NaN;
 
        end
    

        
        
        Hs(cell_lon_idx,wave_lat-i-1) = 4*sqrt(sum(S_attn(i,:).*dwavew));
        Hs_dropoff(cell_lon_idx,wave_lat-i-1) = 4*sqrt(sum(S_attn_dropoff(i,:).*dwavew));

        % Calculate the wave height for each wave type
        %lambda = gravity./(2*pi*(wavefreq).^2); % Wavelength, m
        %idx_wind = lambda < 50;
        %idx_swell = lambda > 50 & lambda < 154; % 154 instead of 150 as one point is at 153
        %idx_long = lambda > 154;
        %Hs_wind(cell_lon_idx,wave_lat-i-1) = 4*sqrt(sum(S_attn_adj(i,idx_wind).*dwavew(idx_wind)));
        %Hs_swell(cell_lon_idx,wave_lat-i-1) = 4*sqrt(sum(S_attn_adj(i,idx_swell).*dwavew(idx_swell)));
        %Hs_long(cell_lon_idx,wave_lat-i-1) = 4*sqrt(sum(S_attn_adj(i,idx_long).*dwavew(idx_long)));
    end % i, prop_length
end


 land_mask = isnan(aice_data(cell_lon_idx,:));

 Hs(cell_lon_idx,land_mask) = NaN;
 Hs_dropoff(cell_lon_idx,land_mask) = NaN;

 

 Hs(Hs==0) = NaN;
 Hs_dropoff(Hs_dropoff==0) = NaN;

 %%
close all 
conFigure(11)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'mlinelimit', [-75 -50]);
    setm(w, 'plinelimit', [-75 -50]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,Hs) 
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    a = colorbar;
    a.Label.String = "$H_s$ [m]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([0,10])
    cmap_cols = cmocean('balance','pivot',0,31);%-10^(-3));
    %cmap_cols = cmap_cols.^2;
    %cmap_temp = cmocean('balance',10,'pivot',0);
    %cmap_cols(end,:) = [0.9,0.9,.9];
    set(gca,'ColorScale','linear')
    colormap(cmap_cols)
    title(filename(end-12:end-3))
    %colormap(cmocean('amp'))
exportgraphics(f,strcat('HsMBK',filename(end-12:end-3),'.pdf'),'ContentType','vector')

conFigure(11)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'mlinelimit', [-75 -50]);
    setm(w, 'plinelimit', [-75 -50]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,Hs_dropoff) 
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    a = colorbar;
    a.Label.String = "$H_s$ [m]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([0,10])
    cmap_cols = cmocean('balance','pivot',0,31);%-10^(-3));
    %cmap_cols = cmap_cols.^2;
    %cmap_temp = cmocean('balance',10,'pivot',0);
    %cmap_cols(end,:) = [0.9,0.9,.9];
    set(gca,'ColorScale','linear')
    colormap(cmap_cols)
    title(filename(end-12:end-3))
    %colormap(cmocean('amp'))
exportgraphics(f,strcat('HsDropoff',filename(end-12:end-3),'.png'),'ContentType','image','Resolution',1080)

conFigure(11)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'mlinelimit', [-75 -50]);
    setm(w, 'plinelimit', [-75 -50]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,Hs_dropoff-Hs) 
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    a = colorbar;
    a.Label.String = "$H_s$ [m]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([-2,2])
    cmap_cols = cmocean('balance','pivot',0,31);%-10^(-3));
    %cmap_cols = cmap_cols.^2;
    %cmap_temp = cmocean('balance',10,'pivot',0);
    %cmap_cols(end,:) = [0.9,0.9,.9];
    set(gca,'ColorScale','linear')
    colormap(cmap_cols)
    title(filename(end-12:end-3))
    %colormap(cmocean('amp'))
exportgraphics(f,strcat('HsDropoff-MBK',filename(end-12:end-3),'.png'),'ContentType','image','Resolution',1080)

%fsdrad_data 
idx = aice_data < 0.01;
fsdrad_data(idx) = NaN;
conFigure(11)
f = figure;
w = worldmap('world');
    axesm eqaazim; %, eqaazim wetch eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
    setm(w, 'Origin', [-90 0 0]);
    setm(w, 'maplatlimit', [-90,-50]);
    setm(w, 'maplonlimit', [-180,180]);
    setm(w, 'meridianlabel', 'off')
    setm(w, 'parallellabel', 'off')
    setm(w, 'mlabellocation', 60);
    setm(w, 'plabellocation', 10);
    setm(w, 'mlabelparallel', -45);
    setm(w, 'mlinelimit', [-75 -50]);
    setm(w, 'plinelimit', [-75 -50]);
    setm(w, 'grid', 'off');
    setm(w, 'frame', 'off');
    setm(w, 'labelrotation', 'on')
    pcolorm(lat,lon,fsdrad_data) 
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
    a = colorbar;
    a.Label.String = "Representative floe radius $r_a$ [m]";
    a.TickLabelInterpreter = 'latex';
    a.Label.Interpreter = 'latex';
    caxis([1,1000])
    cmap_cols = cmocean('haline',31);%-10^(-3));
    %cmap_cols = cmap_cols.^2;
    %cmap_temp = cmocean('balance',10,'pivot',0);
    %cmap_cols(end,:) = [0.9,0.9,.9];
    set(gca,'ColorScale','log')
    colormap(cmap_cols)
    title(filename(end-12:end-3))
    %colormap(cmocean('amp'))
exportgraphics(f,strcat('FSDrad',filename(end-12:end-3),'.png'),'ContentType','image','Resolution',1080)


%%
%% Functions


function [S_out,beta_N] = amplitude_dropoff(S_in,omega,floe_size)
    % Setup
    gravity = 9.80665;
    T = (2*pi)./omega;
    lambda = (gravity*(T.^2))/(2*pi); % Wave lengths, dispersion relation, m
    
    % Scattering calculation
    % Parameters
    max_energy_loss = 1.0;
    min_energy_loss = 0;
    d_gt_lambda = 5;
    unity = 1;
    amplitude = (max_energy_loss-min_energy_loss)/2;

    d_lambda = 2*floe_size./lambda; 
    idx = 2*floe_size < (lambda/2); % Converting from radius to diameter
    beta_N = zeros(1,length(omega));
    beta_N(idx) = 1;
    beta_N(~idx) = 0.5;
    S_out = S_in.*beta_N;
end

function [S,omega,T] = SDF_Bretschneider(Hs,Tm,nw)
    fmin = 1/1000;%1/50;%1/16; % freq min
    fmax = 1;%1/2;%1/6; % freq max
    om1=2*pi*fmin; % ang freqs, rad/s
    om2=2*pi*fmax;
    om_0 = (om2 - om1)/(nw-1); % steps
    gravit = 9.81;
    % Calculate wave numbers and wavelengths
    for lp_i=1:nw
      omega(lp_i)        = om1 + (lp_i-1)*om_0; % Frequency, rad/s
      T(lp_i)            = 2*pi/omega(lp_i); % Period, s
      lam_wtr_in(lp_i)   = gravit*(T(lp_i)^2)/2/pi; % Wave number
      k_wtr_in(lp_i)     = 2*pi/lam_wtr_in(lp_i); % Wave length
    end
    
    tmin = 1/fmax;
    tmax = 1/fmin;
    t1 = tmin;
    t2 = tmax;
    t_0 = (t2-t1)/(nw-1);
    
    
    %clear omega
    % Define steps in terms of period not frequency
    %for lp_i=1:nw
    %   T(lp_i) = t1 + (lp_i-1)*t_0;
    %   omega(lp_i) = 2*pi/T(lp_i);
    %end


    %  omega(lp_i)        = om1 + (lp_i-1)*om_0; % Frequency, rad/s
    %  T(lp_i)            = 2*pi/omega(lp_i); % Period, s
    %  lam_wtr_in(lp_i)   = gravit*(T(lp_i)^2)/2/pi; % Wave number
    %  k_wtr_in(lp_i)     = 2*pi/lam_wtr_in(lp_i); % Wave length
    %end
    



    
    for lp_i = 1:nw
        om_m = 2*pi/Tm; % Peak angular frequency, rad/s
        tau(lp_i)  = 2*pi/omega(lp_i); % s
    end
    
    moment_no = 0;
    f1 = (5/16)*(Hs.^2)*(om_m.^4); % m^2/s^4
    f2 = omega.^(moment_no-5); % rad/s^-5
    f3 = exp(-1.25*((tau/Tm).^4)); % dimensionless (as exponentials have no dimension)

    S = f1.*f2.*f3; % rad x m^2 s/rad
end


function [S,T] = SDF_BretschneiderT(Hs,Tm,nw)
    fmin = 1/400;%1/50;%1/16; % freq min
    fmax = 1;%1/2;%1/6; % freq max

    Tmin = 1/fmax; % s
    Tmax = 1/fmin; % s

    T1=Tmin; % ang freqs, rad/s
    T2=Tmax;
    T_0 = (T2 - T1)/(nw-1); % steps
    gravit = 9.81;
    % Calculate wave numbers and wavelengths
    for lp_i=1:nw
      T(lp_i)        =  T1 + (lp_i-1)*T_0; % Period, s
      %lam_wtr_in(lp_i)   = gravit*(T(lp_i)^2)/2/pi; % Wave length
      %k_wtr_in(lp_i)     = 2*pi/lam_wtr_in(lp_i); % Wave number
    end
    
%    tmin = 1/fmax;
%    tmax = 1/fmin;
%    t1 = tmin;
%    t2 = tmax;
%    t_0 = (t2-t1)/(nw-1);
    
%    for lp_i = 1:nw
%        om_m = 2*pi/Tm; % Peak angular frequency, rad/s
%        tau(lp_i)  = 2*pi/omega(lp_i); % s
%    end
    
    moment_no = 0;
    f1 = (5/16)*(Hs.^2)*(Tm.^(-4)); % m^2/s^4
    f2 = T.^(5); % rad/s^-5
    f3 = exp(-1.25*((T./Tm).^4)); % dimensionless (as exponentials have no dimension)
    S = f1.*f2.*f3; % rad x m^2 s/rad
end

function attn_spec = fn_Attn_MBK(local_wave_spec)
    % Attenuate according to Meylan et al. (2014)
    dum_om = local_wave_spec;
    beta0 = 5.376168295200780e-005;
    beta1 = 2.947870279251530e-005;
    
    fn_Attn_MBK1 = beta0*(dum_om.^2) + beta1*(dum_om.^4);
    attn_fac = 1;
    attn_spec = attn_fac*fn_Attn_MBK1;
end


function [int_E_f_theta,D,theta_vec] = cosine_spreader(S_init,theta_m,thn,n)
    % Spread the wave spectrum through angular space
    % theta 
    % theta_m is the mean wave direction [rad], South = 0
    % and integrate over -pi/2 to pi/2
    % theta0 is the MWD
    % D is the energy in directional spectrum
    % theta_vec is the corresponding angles for D
    % n is the index

    theta_vec = linspace(theta_m-pi,theta_m+pi,thn);
    low_bnd = theta_m - pi < theta_vec;
    upp_bnd = theta_vec < theta_m + pi;
    bnd = low_bnd.*upp_bnd;
    D = zeros(1,thn);
    D = cos((theta_vec-theta_m)/2).^(2*n);
    D(~bnd) = 0; % Limiting the widths of the cosine in [theta_m - pi/2, theta0 + pi/2] 

    dtheta = theta_vec(1:end) - [0,theta_vec(1:end-1)];
    C = 1/sum(D.*dtheta); % Normalising constant
    D = C*D; % Normalising the directional spectrum

    %D = (2/pi).*(cos(theta_vec-theta0)).^n;
    %
    % Integration step
    %upper = pi/2;
    %lower = 3*pi/2;
    %low_bnd = lower < theta_vec;
    %upp_bnd = theta_vec < upper; 
    %%bnd = low_bnd | upp_bnd;
    %D(~bnd) = 0;
    %dtheta = theta_vec(1:end) - [0,theta_vec(1:end-1)];
    %int_D = sum(D.*dtheta);
    
    % Integration step
%    low_bnd = theta0 - pi/2 < theta_vec;
%    upp_bnd = theta_vec < theta0 + pi/2;
%    bnd = low_bnd.*upp_bnd;
%    D(~bnd) = 0;
%    dtheta = theta_vec(1:end) - [0,theta_vec(1:end-1)];
%    int_D = sum(D.*dtheta);
%    D = D./int_D;

    % Integrate over the Southern wedge
    lower = -pi/2;
    upper = pi/2;
    low_bnd = lower < theta_vec;
    upp_bnd = theta_vec < upper;
    bnd = low_bnd & upp_bnd;
    lower = 3*pi/2;
    upper = pi/2;
    low_bnd = lower < theta_vec;
    upp_bnd = theta_vec < upper;
    bnd2 = low_bnd | upp_bnd;
    bnd = bnd | bnd2;
    E_f_theta = S_init'*D;
    int_E_f_theta = sum(S_init'*(dtheta(bnd).*D(bnd)),2)';

end

function [data_out] = data_format_cawcr(filedir,variable)
    % OM2 grid

    lon = ncread(filedir,'LON');
    lat = ncread(filedir,'LAT');

    row = 281;
    dim = 3;
    data1 = ncread(filedir, variable);
    
    latitude = [-90,90];
    longitude = [-180,180];
    

    data1 = rearrange_matrix(data1,row,dim);
    
    % fixing data
    [m, ~] = size(lon);
    lon = [lon; lon(end,:) + 360/m];
    lat = [lat; lat(end,:)];
    data_out = [data1; data1(end,:,:)];

end




function [S_out] = wave_attenuation(conc,L,nw,S_in,omega,option,floe_size)
% conc is aice
% L is the length of the cell [m]
% nw is the number of points in frequency space
% S_in is the incoming wave spectra
% option : integer
     % Do MBK exponential attenuation
    if option == "MBK"
    for lp_i=1:nw
      alpha(lp_i)  = conc*fn_Attn_MBK(omega(lp_i))/0.7; % 0.7 comes from the concentration measures in the observations of MBK
      S_out(lp_i) = S_in(lp_i)*exp(-alpha(lp_i)*L);
    end 
    else
        % Adjusted alpha
        gravity = 9.80665;
        T = (2*pi)./omega;
        lambda = (gravity*(T.^2))/(2*pi); % Wave lengths
    
        % Parameters
        max_alpha = 2.0;
        min_alpha = 0.5;
        unity = 1;
        d_gt_lambda = 2;

        d_lambda = 2*floe_size./lambda; % Diameter/wavelength
        vert_translation = (max_alpha+min_alpha)/2;
        horz_translation = (d_gt_lambda+unity)/2;
        amplitude = (max_alpha-min_alpha)/2;
        
        for lp_i=1:nw
            alpha_M(lp_i) = conc*fn_Attn_MBK(omega(lp_i))/0.7; % 0.7 comes from the concentration measures in the observations of MBK
            alpha_N(lp_i) = amplitude*tanh(-(d_lambda(lp_i)-horz_translation)) + vert_translation;
            alpha(lp_i)  = alpha_N(lp_i)*alpha_M(lp_i); 
            S_out(lp_i) = S_in(lp_i)*exp(-alpha(lp_i)*L);
        end 
    end
end
