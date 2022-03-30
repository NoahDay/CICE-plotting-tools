%% FSTD_VALIDATION
% This script will validate the FSTD construction through CICE history
% fields to make sure everything is behaving as expected.
% Aim: To get from afsdn to aice
% Method:
% 1. Integrate aicen to aice
% 2. Integrate 


clear all
close all
addpath functions
addpath packages/bedmap2_toolbox_v4
filename = 'cases/ocntest/history/iceh.2009-01-01.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';
ssd = 0;
sector = "SH";

%% Preamble
close all
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'wimoninit';
sector = "SH";
ssd = 1;
if ssd == 1
    ssd_dir = '/Volumes/NoahDay5TB/cases/';
    filedir = strcat(ssd_dir,case_name);
else
    filedir = strcat('cases/',case_name);
end
grid = 'gx1'; 

day = 31;
month_init = 7;
year = 2005;
date = sprintf('%d-0%d-%d', year, month_init, day);

figcount = 0;
filename = strcat(filedir,"/history/iceh.",date,".nc");

[lat, lon, row] = grid_read(grid);

NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;
thick_rad_l = [0; NCAT(1:end-1)]; % Thickness lower bound
thick_rad_h = NCAT(1:end); % Floe radius higher bound
thick_binwidth = thick_rad_h - thick_rad_l; 
thick_binwidth(end) = 15;


% 1. Integrate aicen wrt to ice thickness to get aice
aicen_data = data_format_sector(filename,"aicen",sector);
aice_data = data_format_sector(filename,"aice",sector);

[nx,ny,~] = size(aicen_data);
% a = INT(g(h)dh)
for i = 1:nx
    for j = 1:ny
        temp(:) = aicen_data(i,j,:);
        aice_transformed(i,j) = sum(temp'.*thick_binwidth);
        clear temp
    end
end
t1 = tiledlayout(1,3);
nexttile
map_plot(aice_data,"aice",sector,grid,[0,1]);
title("Model data")
nexttile
map_plot(aice_transformed,"aice",sector,grid,[0,1]);
title("Transformed data - fix bin width")
nexttile
map_plot(aice_data-aice_transformed,"aice",sector,grid,[-1,1]);
title("Difference")


%% 2. Integrate afsd wrt to floe size to get aice
clear aice_transformed
close all
afsd_data = data_format_sector(filename,"afsd",sector);
aice_data = data_format_sector(filename,"aice",sector);

[nx,ny,~] = size(afsd_data);
% a = INT(g(h)dh)
for i = 1:nx
    for j = 1:ny
        temp(:) = afsd_data(i,j,:);
        if isnan(temp)
            aice_transformed(i,j) = NaN;
        else
            aice_transformed(i,j) = sum(temp.*floe_binwidth);
        end
        clear temp
    end
end
f2 = figure(2);
t2 = tiledlayout(1,3);
t2.TileSpacing = 'compact';
fontsize = 14;


nexttile
[p,a] = map_plot(aice_data,"aice",sector,grid,[0,1]);
t = title("\textbf{CICE output}");
t.Interpreter = "latex";
t.FontSize = fontsize;
a.FontSize = fontsize;
a.TickLabelInterpreter = "Latex";
cmocean('ice') 

nexttile
[p,a] = map_plot(aice_transformed,"aice",sector,grid,[0,1]);
t = title("\textbf{$\sum_{j=1}^{n_f} F(r_j)\Delta r_j$}");
t.Interpreter = "latex";
t.FontSize = fontsize;
a.FontSize = fontsize;
a.TickLabelInterpreter = "latex";
cmocean('ice') 

nexttile
[p,a] = map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.05,0.05]);
cmocean('balance') 
t = title("\textbf{Difference}");
t.Interpreter = "latex";
t.FontSize = fontsize;
a.Label.String = 'Ice concentration';
a.FontSize = fontsize;
a.TickLabelInterpreter = "latex";
f2.Position = [800 1000 1000 400];

%exportgraphics(f2,'aice.pdf','ContentType','vector')

icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);
aice_transformed_vec = aice_transformed(icemask);
error_vec = aice_vec-aice_transformed_vec;
max(abs(error_vec))


f3 = figure(3);
t = tiledlayout(2,2);
t.TileSpacing = 'compact';


nexttile
% Correlation plot
scatter(aice_transformed_vec, aice_vec)
hold off
ax = gca;
ax.FontSize = fontsize; 
y = ylabel('CICE output $a_i$');
y.FontSize = fontsize;
x = xlabel('My calculation of $a_i$');
x.FontSize = fontsize;


nexttile
% Box plot comparison for error and actual value
set(groot, 'defaultAxesTickLabelInterpreter','latex');
 boxplot([error_vec,aice_vec],'labels',{'Error','$a_i$'});
ax = gca;
ax.FontSize = fontsize; 
ax.TickLabelInterpreter='latex';
ylabel("Ice concentration",'Interpreter','latex')

nexttile
scatter(aice_vec,error_vec)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); %trying to set the default
ax = gca;
ax.FontSize = fontsize; 
xlabel("Ice concentration")
ylabel("Error (concentration)")

nexttile
scatter(lat(icemask),error_vec)
ax = gca;
ax.FontSize = fontsize; 
xlabel("Latitude",'Interpreter','latex')
ylabel("Error (concentration)",'Interpreter','latex')

%exportgraphics(f3,'aice_error.pdf','ContentType','vector')

% Statistics
fprintf("Data statistics: \n")
fprintf("Correlation between methods: %g\n", corr(aice_transformed_vec, aice_vec))
fprintf('The max fsdrad is: %d\n',max(aice_vec))
fprintf('The mean fsdrad is: %d\n',mean(aice_vec))
fprintf('The median fsdrad is: %d\n',median(aice_vec))
fprintf('The mode fsdrad is: %d\n',mode(aice_vec))
fprintf('The std fsdrad is: %d\n',std(aice_vec))
fprintf("Error statistics: \n")
fprintf('The max error is: %d\n',max(error_vec))
fprintf('The mean error is: %d\n',mean(error_vec))
fprintf('The median error is: %d\n',median(error_vec))
fprintf('The mode error is: %d\n',mode(fsdrad_vec))
fprintf('The std error is: %d\n',std(error_vec))






% 3. Integrate afsdn to get aicen


%% 4. Integrate afsdn to get afsd
clear afsd_transformed 
close all
afsd_data(:,:,:) = data_format_sector(filename,"afsd",sector);
afsd_transformed(:,:,:) = fsd_converter(filename,"afsdn","afsd");

figcount = figcount + 1;
f4 = figure(figcount);
t4 = tiledlayout(2,8);
for i = 1:Nf
    nexttile
    lims = max(max(afsd_data(:,:,i) - afsd_transformed(:,:,i)));
    [p,a] = map_plot(afsd_data(:,:,i) - afsd_transformed(:,:,i),"aice",sector,grid,[-lims,lims]);
    t = title(sprintf("FSD Cat %d",i));
    cmocean('balance') 
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    %a.Label.String = 'AFSD';
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "latex";
end
t= title(t4,"F(r) - $\sum_{i=1}^{n_c} f(r,h) $");
t.Interpreter = "latex";
t.FontSize = fontsize;
f4.Position = [800 1000 1200 400];
%exportgraphics(f4,'afsd.pdf','ContentType','vector')
%%

icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);

figcount = figcount + 1;
f5 = figure(figcount);
t5 = tiledlayout(2,8);
for i = 1:Nf
    nexttile
    temp_trans = afsd_transformed(:,:,i);
    afsd_transformed_vec = temp_trans(icemask);
    temp_raw = afsd_data(:,:,i);
    afsd_raw_vec = temp_raw(icemask);
    
    error_vec = abs(afsd_raw_vec-afsd_transformed_vec);
    scatter(aice_vec,error_vec)
    title(sprintf("FSD Cat %d",i))
    xlabel("Ice concentration")
    ylabel("Error")
    clear temp_trans temp_raw
end

figcount = figcount + 1;
f6 = figure(figcount);
t6 = tiledlayout(2,8);
for i = 1:Nf
    nexttile
    temp_trans = afsd_transformed(:,:,i);
    afsd_transformed_vec = temp_trans(icemask);
    temp_raw = afsd_data(:,:,i);
    afsd_raw_vec = temp_raw(icemask);
    
    error_vec = abs(afsd_raw_vec-afsd_transformed_vec);
    scatter(lat(icemask),error_vec)
    title(sprintf("FSD Cat %d",i))
    clear temp_trans temp_raw
end

% Box plot comparison for error and actual value
figcount = figcount + 1;
f7 = figure(figcount);
t7 = tiledlayout(2,8);
for i = 1:Nf
    nexttile
    temp_trans = afsd_transformed(:,:,i);
    afsd_transformed_vec = temp_trans(icemask);
    temp_raw = afsd_data(:,:,i);
    afsd_raw_vec = temp_raw(icemask);
    
    error_vec = abs(afsd_raw_vec-afsd_transformed_vec);
    boxplot([error_vec,afsd_raw_vec])
    title(sprintf("FSD Cat %d",i))
    if i == 1
        %legend({"Error","AFSD"})
    end
    clear temp_trans temp_raw
end


%% 5. Get FSDrad from AFSDN
clear fsdrad_transformed 
clc
close all
fsdrad_data(:,:) = data_format_sector(filename,"fsdrad",sector);
fsdrad_transformed(:,:) = fsd_converter(filename,"afsdn","fsdrad");

f2 = figure(2);
t2 = tiledlayout(1,3);
nexttile
[p,a] = map_plot(fsdrad_data,"fsdrad",sector);
t = title("$r_a$");
t.Interpreter = "latex";
t.FontSize = fontsize;
a.FontSize = fontsize;
a.TickLabelInterpreter = "latex";
cmocean('ice') 
nexttile
[p,a] = map_plot(fsdrad_transformed,"fsdrad",sector);
t = title("$r_a$ calculated from \texttt{afsdn}");
t.Interpreter = "latex";
t.FontSize = fontsize;

a.FontSize = fontsize;
a.TickLabelInterpreter = "latex";
cmocean('ice') 
nexttile
lims = max(max(fsdrad_data-fsdrad_transformed));
[p,a] = map_plot(fsdrad_data-fsdrad_transformed,"fsdrad",sector,grid,[-lims,lims]);
t = title("Difference");
t.Interpreter = "latex";
t.FontSize = fontsize;
a.FontSize = fontsize;
a.Label.String = 'representative radius (m)';
a.TickLabelInterpreter = "latex";
cmocean('balance') 
f2.Position = [800 1000 1200 400];
%exportgraphics(f2,'fsdrad.pdf','ContentType','vector')
%%
icemask = aice_data > 0.01;
fsdrad_vec = fsdrad_data(icemask);
fsdrad_transformed_vec = fsdrad_transformed(icemask);
error_vec = fsdrad_vec-fsdrad_transformed_vec;

f1 = figure(1);
t = tiledlayout(2,2);
t.TileSpacing = 'compact';


nexttile
% Correlation plot
scatter(fsdrad_transformed_vec, fsdrad_vec)
hold off
ax = gca;
ax.FontSize = fontsize; 
y = ylabel('CICE output $r_a$ (m)');
y.FontSize = fontsize;
x = xlabel('My calculation of $r_a$ (m)');
x.FontSize = fontsize;


nexttile
% Box plot comparison for error and actual value
set(groot, 'defaultAxesTickLabelInterpreter','latex');
 boxplot([error_vec,fsdrad_vec],'labels',{'Error','$r_a$'});
ax = gca;
ax.FontSize = fontsize; 
ax.TickLabelInterpreter='latex';
ylabel("Floe size (m)",'Interpreter','latex')

nexttile
scatter(aice_vec,error_vec)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); %trying to set the default
ax = gca;
ax.FontSize = fontsize; 
xlabel("Ice concentration")
ylabel("Error (m)")

nexttile
scatter(lat(icemask),error_vec)
ax = gca;
ax.FontSize = fontsize; 
xlabel("Latitude",'Interpreter','latex')
ylabel("Error (m)",'Interpreter','latex')

%exportgraphics(f1,'error_fsdrad.pdf','ContentType','vector')



% Statistics
fprintf("Data statistics: \n")
fprintf("Correlation between methods: %g\n", corr(fsdrad_transformed_vec, fsdrad_vec))
fprintf('The max fsdrad is: %d\n',max(fsdrad_vec))
fprintf('The mean fsdrad is: %d\n',mean(fsdrad_vec))
fprintf('The median fsdrad is: %d\n',median(fsdrad_vec))
fprintf('The mode fsdrad is: %d\n',mode(fsdrad_vec))
fprintf('The std fsdrad is: %d\n',std(fsdrad_vec))
fprintf("Error statistics: \n")
fprintf('The max error is: %d\n',max(error_vec))
fprintf('The mean error is: %d\n',mean(error_vec))
fprintf('The median error is: %d\n',median(error_vec))
fprintf('The mode error is: %d\n',mode(fsdrad_vec))
fprintf('The std error is: %d\n',std(error_vec))



%% Functions
function [p] = plot_error(fsdrad_transformed_vec,fsdrad_vec)
    t = tiledlayout(2,2);
    t.TileSpacing = 'compact';

    nexttile
    % Correlation plot
    %f6 = figure(6);
    scatter(fsdrad_transformed_vec, fsdrad_vec)
    hold off
    %l = legend({'Error','$r_a$'});
    %l.FontSize = fontsize;
    ax = gca;
    ax.FontSize = fontsize; 
    y = ylabel('CICE output $r_a$ (m)');
    y.FontSize = fontsize;
    x = xlabel('My calculation of $r_a$ (m)');
    x.FontSize = fontsize;
    %exportgraphics(f6,'histfsdrad.pdf','ContentType','vector')

    error_vec = fsdrad_vec - fsdrad_transformed_vec;
    
    nexttile
    % Box plot comparison for error and actual value
    %f5 = figure(5);
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
     boxplot([error_vec,fsdrad_vec],'labels',{'Error','$r_a$'});
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    ylabel("Floe size (m)",'Interpreter','latex')
    %exportgraphics(f5,'boxplotfsdrad.pdf','ContentType','vector')

    nexttile
    %f3 = figure(3);
    scatter(aice_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Ice concentration")
    ylabel("Error (m)")
    %exportgraphics(f3,'err_conc.pdf','ContentType','vector')

    nexttile
    %f4 = figure(4);
    scatter(lat(icemask),error_vec)
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Latitude",'Interpreter','latex')
    ylabel("Error (m)",'Interpreter','latex')

end

