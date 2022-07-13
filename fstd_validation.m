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
filename = '/Users/noahday/GitHub/CICE-plotting-tools/cases/pancake_tracer/history/iceh.2005-09.nc';
% Read the header
ncdisp(filename)
% 
grid = 'gx1';
sector = "SH";

%% Preamble
close all
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'forcingoff';
sector = "SH";
ssd = 0;
if ssd == 1
    ssd_dir = '/Volumes/NoahDay5TB/cases/';
    filedir = strcat(ssd_dir,case_name);
else
    filedir = strcat('cases/',case_name);
end
grid = 'gx1'; 

day = 10;
month_init = 9;
year = 2005;
date = sprintf('%d-0%d-%d', year, month_init, day);

figcount = 0;
%filename = strcat(filedir,"/history/iceh.",date,".nc");

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


%% 1. Integrate aicen wrt to ice thickness to get aice
clear temp
aicen_data = data_format_sector(filename,"aicen",sector);
vicen_data = data_format_sector(filename,"vicen",sector);
aice_data = data_format_sector(filename,"aice",sector);

[nx,ny,~] = size(aicen_data);
% a = INT(g(h)dh)
for i = 1:nx
    for j = 1:ny
        temp(:) = aicen_data(i,j,:);
        thick_binwidth = 1;
        aice_transformed(i,j) = sum(temp');
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
map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.1,0.1]);
title("Difference")


icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);
aice_transformed_vec = aice_transformed(icemask);
error_vec = aice_vec-aice_transformed_vec;
max(abs(error_vec))

% Statistics
fprintf("Data statistics: \n")
fprintf("Correlation between methods: %g\n", corr(aice_transformed_vec, aice_vec))
fprintf('The max aice is: %d\n',max(aice_vec))
fprintf('The mean aice is: %d\n',mean(aice_vec))
fprintf('The median aice is: %d\n',median(aice_vec))
fprintf('The mode aice is: %d\n',mode(aice_vec))
fprintf('The std aice is: %d\n',std(aice_vec))
fprintf("Error statistics: \n")
fprintf('The max error is: %d\n',max(error_vec))
fprintf('The mean error is: %d\n',mean(error_vec))
fprintf('The median error is: %d\n',median(error_vec))
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))



%% 2. Integrate afsd wrt to floe size to get aice
clear aice_transformed f1 f2
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


% Plot comparison
f1 = figure(1);
t1 = tiledlayout(1,3);
t1.TileSpacing = 'compact';
fontsize = 14;

nexttile
    [p,a] = map_plot(aice_data,"aice",sector,grid,[0,1]);
    t = title(" \textbf{\texttt{aice} CICE output}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "Latex";
    cmocean('ice',10) 

nexttile
    [p,a] = map_plot(aice_transformed,"aice",sector,grid,[0,1]);
    t = title("\textbf{\texttt{aice} calculated from \texttt{afsd}}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "latex";
    cmocean('ice',10) 
    clear a

nexttile
    [p,a] = map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.05,0.05]);
    cmocean('balance',15) 
    t = title("\textbf{Difference}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.Label.String = 'Ice concentration';
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "Latex";
    f1.Position = [800 1000 1000 400];

%exportgraphics(f1,'aice.pdf','ContentType','vector')

icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);
aice_transformed_vec = aice_transformed(icemask);
error_vec = aice_vec-aice_transformed_vec;
max(abs(error_vec))

% Plot errors
f2 = figure(2);
t = tiledlayout(2,2);
t.TileSpacing = 'compact';


nexttile
    % Correlation plot
    scatter(aice_transformed_vec, aice_vec)
    hold on
    plot([0,1],[0,1])
    hold off
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    y = ylabel('CICE output $a_i$','Interpreter','latex');
    y.FontSize = fontsize;
    x = xlabel('My calculation of $a_i$','Interpreter','latex');
    x.FontSize = fontsize;


nexttile
    % Histogram comparison
    [h1,edges] = histcounts(aice_vec, 10);
    [h2,edges] = histcounts(aice_transformed_vec, 10);
    ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
    bar(ctrs, log([h1 ;h2])')
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Ice concentration",'Interpreter','latex')
    ylabel("log(Counts)",'Interpreter','latex')
    legend("CICE","Calculated",'Location','northwest')
    

nexttile
    % Error vs aice
    scatter(aice_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    yline(0)
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Ice concentration")
    ylabel("Error (concentration)")

nexttile
   % Box plot comparison for error and actual value
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
     boxplot([error_vec,aice_vec],'labels',{'Error','$a_i$'});
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    ylabel("Ice concentration",'Interpreter','latex')

%exportgraphics(f2,'aice_error.pdf','ContentType','vector')

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
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))



%% Make latex table
%variables_names = ['Mean';'Median';'Mode';'Standard deviation';'Maximum'];
Mean = [mean(aice_vec), mean(error_vec)];
%patients = table(LastName,Age,Smoker,Height,Weight,BloodPressure)
data_table = table(Mean);
err = stat_vec(error_vec);
%tex_data = [variables_names,vpa(stat_vec(error_vec))];
%latex([variables_names,vpa(stat_vec(error_vec))])

latex(['Mean',vpa(round(mean(aice_vec),5)), vpa(round(mean(error_vec),5));...
       'Median',vpa(round(median(aice_vec),5)), vpa(round(median(error_vec),5));...
       'Mode',vpa(round(mode(aice_vec),5)), vpa(round(mode(error_vec),5));...
       'STD',vpa(round(std(aice_vec),5)), vpa(round(std(error_vec),5));...
       'Maximum',vpa(round(max(aice_vec),5)), vpa(round(max(error_vec),5));...
    ])


%% 3. Integrate afsdn to get aicen 
close all
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
afsdn_data(:,:,:,:) = data_format_sector(filename,"afsdn",sector);
for i = 1:Nc
    for j = 1:nx
        for k = 1:ny
            temp(:) = afsdn_data(j,k,:,i);
            aicen_transformed(j,k,i) =  sum(temp.*floe_binwidth);
        end
    end
end

f1 = figure;
t1 = tiledlayout(1,5);
for i = 1:Nc
    nexttile
    lims = max(max(aicen_data(:,:,i) - aicen_transformed(:,:,i)));
    [p,a] = map_plot(aicen_data(:,:,i) - aicen_transformed(:,:,i),"aice",sector,grid,[-lims,lims]);
    t = title(sprintf("ITD Cat %d",i));
    cmocean('balance') 
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    %a.Label.String = 'AFSD';
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "latex";
end
t= title(t1,"g(h) - $\sum_{i=1}^{n_f} f(r_i,h) $");
t.Interpreter = "latex";
t.FontSize = fontsize;
f1.Position = [800 1000 1200 400];
%exportgraphics(f4,'afsd.pdf','ContentType','vector')


% Integrate again to get aice
aice_data = data_format_sector(filename,"aice",sector);

[nx,ny,~] = size(aicen_data);
clear temp
% a = INT(g(h)dh)
for i = 1:nx
    for j = 1:ny
        temp(:) = aicen_transformed(i,j,:);
        thick_binwidth = 1;
        aice_transformed(i,j) = sum(temp');
        clear temp
    end
end
f2 = figure;
t2 = tiledlayout(1,3);
nexttile
map_plot(aice_data,"aice",sector,grid,[0,1]);
title("Model data")
nexttile
map_plot(aice_transformed,"aice",sector,grid,[0,1]);
title("Transformed data - fix bin width")
nexttile
map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.1,0.1]);
title("Difference")


icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);
aice_transformed_vec = aice_transformed(icemask);
error_vec = aice_vec-aice_transformed_vec;
max(abs(error_vec))

% Statistics
fprintf("Data statistics: \n")
fprintf("Correlation between methods: %g\n", corr(aice_transformed_vec, aice_vec))
fprintf('The max aice is: %d\n',max(aice_vec))
fprintf('The mean aice is: %d\n',mean(aice_vec))
fprintf('The median aice is: %d\n',median(aice_vec))
fprintf('The mode aice is: %d\n',mode(aice_vec))
fprintf('The std aice is: %d\n',std(aice_vec))
fprintf("Error statistics: \n")
fprintf('The max error is: %d\n',max(error_vec))
fprintf('The mean error is: %d\n',mean(error_vec))
fprintf('The median error is: %d\n',median(error_vec))
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))



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

%% and then integrate "afsd" to get aice

clear aice_transformed f1 f2
close all
aice_transformed(:,:) = fsd_converter(filename,"afsdn","aice");
aice_data = data_format_sector(filename,"aice",sector);


% Plot comparison
f1 = figure(1);
t1 = tiledlayout(1,3);
t1.TileSpacing = 'compact';
fontsize = 14;

nexttile
    [p,a] = map_plot(aice_data,"aice",sector,grid,[0,1]);
    t = title(" \textbf{\texttt{aice} CICE output}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "Latex";
    cmocean('ice',10) 

nexttile
    [p,a] = map_plot(aice_transformed,"aice",sector,grid,[0,1]);
    t = title("\textbf{\texttt{aice} calculated from \texttt{afsdn}}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "latex";
    cmocean('ice',10) 
    clear a

nexttile
    [p,a] = map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.05,0.05]);
    cmocean('balance',15) 
    t = title("\textbf{Difference}");
    t.Interpreter = "latex";
    t.FontSize = fontsize;
    a.Label.String = 'Ice concentration';
    a.FontSize = fontsize;
    a.TickLabelInterpreter = "Latex";
    f1.Position = [800 1000 1000 400];

%exportgraphics(f1,'afsdn_aice.pdf','ContentType','vector')

icemask = aice_data > 0.01;
aice_vec = aice_data(icemask);
aice_transformed_vec = aice_transformed(icemask);
error_vec = aice_vec-aice_transformed_vec;
max(abs(error_vec))

% Plot errors
f2 = figure(2);
t = tiledlayout(2,2);
t.TileSpacing = 'compact';


nexttile
    % Correlation plot
    scatter(aice_transformed_vec, aice_vec)
    hold on
    plot([0,1],[0,1])
    hold off
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    y = ylabel('CICE output $a_i$');
    y.FontSize = fontsize;
    x = xlabel('My calculation of $a_i$');
    x.FontSize = fontsize;


nexttile
    % Histogram comparison
    [h1,edges] = histcounts(aice_vec, 10);
    [h2,edges] = histcounts(aice_transformed_vec, 10);
    ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
    bar(ctrs, log([h1 ;h2])')
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Ice concentration",'Interpreter','latex')
    ylabel("log(Counts)",'Interpreter','latex')
    legend("CICE","Calculated",'Location','northwest')
    

nexttile
    % Error vs aice
    scatter(aice_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    yline(0)
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Ice concentration")
    ylabel("Error (concentration)")

nexttile
   % Box plot comparison for error and actual value
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
     boxplot([error_vec,aice_vec],'labels',{'Error','$a_i$'});
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    ylabel("Ice concentration",'Interpreter','latex')

exportgraphics(f2,'afsdn_aice_error.pdf','ContentType','vector')

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
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))



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

icemask = aice_data > 0.01;
fsdrad_vec = fsdrad_data(icemask);
fsdrad_transformed_vec = fsdrad_transformed(icemask);
error_vec = fsdrad_vec-fsdrad_transformed_vec;

f1 = figure(1);
t = tiledlayout(2,3);
t.TileSpacing = 'compact';


nexttile
    % Correlation plot
    scatter(fsdrad_transformed_vec, fsdrad_vec)
    hold on
    plot([0,1],[0,1])
    hold off
    ax = gca;
    ax.FontSize = fontsize; 
    y = ylabel('CICE output $r_a$ (m)');
    y.FontSize = fontsize;
    x = xlabel('My calculation of $r_a$ (m)');
    x.FontSize = fontsize;

nexttile
 % Histogram comparison
    [h1,edges] = histcounts(fsdrad_vec, 10);
    [h2,edges] = histcounts(fsdrad_transformed_vec, 10);
    ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
    bar(ctrs, log([h1 ;h2])')
    ax = gca;
    ax.FontSize = fontsize-2; 
    xtickangle(ax,45)
    xlabel("Floe size radius (m)",'Interpreter','latex')
    ylabel("log(Counts)",'Interpreter','latex')
    legend("CICE","Calculated",'Location','north')
    
nexttile([2 1]);
    % Box plot comparison for error and actual value
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
     boxplot([error_vec,fsdrad_vec],'labels',{'Error','$r_a$'});
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    ylabel("Floe size (m)",'Interpreter','latex')      
    
nexttile
    scatter(fsdrad_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    yline(0)
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Floe size (m)")
    ylabel("Error (m)")
    
nexttile
    scatter(aice_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    ax = gca;
    ax.FontSize = fontsize; 
    yline(0)
    xlabel("Ice concentration")
    ylabel("Error (m)")
f1.Position = [800 1000 1200 400];
  

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
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))

%% 5. a) Get FSDrad from AFSD
clear fsdrad_transformed 
clc
close all
fsdrad_data(:,:) = data_format_sector(filename,"fsdrad",sector);
fsdrad_transformed(:,:) = fsd_converter(filename,"afsd","fsdrad");

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
t = title("$r_a$ calculated from \texttt{afsd}");
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

icemask = aice_data > 0.01;
fsdrad_vec = fsdrad_data(icemask);
fsdrad_transformed_vec = fsdrad_transformed(icemask);
error_vec = fsdrad_vec-fsdrad_transformed_vec;

f1 = figure(1);
t = tiledlayout(2,3);
t.TileSpacing = 'compact';


nexttile
    % Correlation plot
    scatter(fsdrad_transformed_vec, fsdrad_vec)
    hold on
    plot([0,1],[0,1])
    hold off
    ax = gca;
    ax.FontSize = fontsize; 
    y = ylabel('CICE output $r_a$ (m)');
    y.FontSize = fontsize;
    x = xlabel('My calculation of $r_a$ (m)');
    x.FontSize = fontsize;

nexttile
 % Histogram comparison
    [h1,edges] = histcounts(fsdrad_vec, 10);
    [h2,edges] = histcounts(fsdrad_transformed_vec, 10);
    ctrs = edges(1)+(1:length(edges)-1).*diff(edges);   % Create Centres
    bar(ctrs, log([h1 ;h2])')
    ax = gca;
    ax.FontSize = fontsize-2; 
    xtickangle(ax,45)
    xlabel("Floe size radius (m)",'Interpreter','latex')
    ylabel("log(Counts)",'Interpreter','latex')
    legend("CICE","Calculated",'Location','north')
    
nexttile([2 1]);
    % Box plot comparison for error and actual value
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
     boxplot([error_vec,fsdrad_vec],'labels',{'Error','$r_a$'});
    ax = gca;
    ax.FontSize = fontsize; 
    ax.TickLabelInterpreter='latex';
    ylabel("Floe size (m)",'Interpreter','latex')      
    
nexttile
    scatter(fsdrad_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    yline(0)
    ax = gca;
    ax.FontSize = fontsize; 
    xlabel("Floe size (m)")
    ylabel("Error (m)")
    
nexttile
    scatter(aice_vec,error_vec)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    ax = gca;
    ax.FontSize = fontsize; 
    yline(0)
    xlabel("Ice concentration")
    ylabel("Error (m)")
f1.Position = [800 1000 1200 400];
  

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
fprintf('The mode error is: %d\n',mode(error_vec))
fprintf('The std error is: %d\n',std(error_vec))


%% 6. Intergrate FSTD wrt fsd to get aicen
clear afsd_transformed 
close all
afsdn_data(:,:,:,:) = data_format_sector(filename,"afsdn",sector);

[nx,ny,~] = size(afsd_data);
% a = INT(g(h)dh)
for i = 1:nx
    for j = 1:ny
        for n = 1:Nc
            temp(:) = afsdn_data(i,j,:,n);
            if isnan(temp)
                aicen_transformed(i,j,n) = NaN;
            else
                aicen_transformed(i,j,n) = sum(temp.*floe_binwidth);
            end
            clear temp
        end
        temp(:) = aicen_transformed(i,j,:);
        thick_binwidth = 1;
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
title("aice from afsdn via aicen")
nexttile
map_plot(aice_data-aice_transformed,"aice",sector,grid,[-0.1,0.1]);
title("Difference")

%% 7. Number FSD
% Get the NFSD from AFSDN
%close all
clear xticks yticks f5 icemask nfsd trcrn
aice_data = data_format_sector(filename,"aice",sector);
icemask0 = aice_data > -eps;
icemask = aice_data > 0.01;
icemask90 = aice_data > 0.9;
%figure(1)
%map_plot(1.0*icemask,"aice",sector,grid,[0,1]);
%figure(2)
%map_plot(1.0*icemask90,"aice",sector,grid,[0,1]);
afsdn_data(:,:,:,:) = data_format_sector(filename,"afsdn",sector);
afsd_data(:,:,:) = data_format_sector(filename,"afsd",sector);
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
tarea = data_format(filename,"tarea");

[nx,ny,nf,nc] = size(afsdn_data);
%icemask = icemask0;


% Make a mask for the Southern Hemisphere
% Coordinates
ocean_mask = data_format(filename,'tmask');
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
sector_mask = false(nx,ny);
for i = 0:lat_out(1)-lat_out(2) % Cycle through latitudes
    sector_mask(lon_out(1):lon_out(3), i+1) = true;
end
sector_mask = logical(sector_mask.*ocean_mask);

% Get tracer array
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            for n = 1:nc
                if aicen_data(i,j,n) < eps
                    trcrn(i,j,k,n) = 0;
                    trcrn2(i,j,k,n) = 0;
                else
                    trcrn(i,j,k,n) = afsdn_data(i,j,k,n)*floe_binwidth(k)*(10^(-3))/aicen_data(i,j,n); % something x Length: (something m) now oonverted to km
                    trcrn2(i,j,k,n) = afsdn_data(i,j,k,n)*floe_binwidth(k)*(10^(-3))/aice_data(i,j);
                end
            end
            if aice_data(i,j) < eps
                afsd_temp(i,j,k) = 0;
            else
                temp(i,j,k) =  sum(afsdn_data(i,j,k,:));
                afsd_temp(i,j,k) =temp(i,j,k)*floe_binwidth(k)*(10^(-3))/aice_data(i,j);
            end
        end
    end
end

% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*(floe_rad_c*(10^(-3))).^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            for n = 1:nc
                % Get the number of floes per cell
                nfstd(i,j,k,n) = trcrn(i,j,k,n)/(floe_area_c(k)); % Dimensionless: something / Length
                nfstd2(i,j,k,n) = trcrn2(i,j,k,n)/(floe_area_c(k)); % Dimensionless: something / Length
                
            end
            nfstd3(i,j,k) = afsd_temp(i,j,k)/(floe_area_c(k)); % Dimensionless: something / Length
        end
    end
end

% Intergrate nfstd to get nfsd

for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            nfsd(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
            nfsd2(i,j,k) = sum(nfstd2(i,j,k,:)); % number of floes per m^2
            nfsd3(i,j,k) = nfstd3(i,j,k); % number of floes per m^2
        end
         nfsd_cell(i,j) = sum(nfsd(i,j,:));
    end
end

for k = 1:nf
    temp = nfsd(:,:,k);
    nfsd_vec(k) = mean(temp(sector_mask)); % Take the cells over the sector
    temp = nfsd2(:,:,k);
    nfsd_vecai(k) = mean(temp(sector_mask)); % Take the cells over the sector
    temp = nfsd3(:,:,k);
    nfsd_vecint(k) = mean(temp(sector_mask)); % Take the cells over the sector
end
%
num_dafsd = numberfsdconverter(filename,afsd_data);

roach_data_sep = [0.3*10^4, 3*10^1, 10^0, 0.8*10^(-1), 0.9*10^(-2), 1.2*10^(-3), 0.5*10^(-4), 0.7*10^(-5), 10^(-6), 1.3*10^(-7), 10^(-7), 10^(-4)];

f1 = figure(1);

t5 = tiledlayout(1,1);%(5,2);
t5.TileSpacing = 'compact';
cust_bounds =  max(NFSD);
xtick = 10.^(0:6);
xticklab = cellstr(num2str(round(log10(xtick(:))), '$10^{%d}$'));
ytick = 10.^(-9:2:13);
yticklab = cellstr(num2str(round(log10(ytick(:))), '$10^{%d}$'));

nfsd_vec2 = nfsd_vec(1:12-1);
nfsd_vec2(12) = sum(nfsd_vec(12:end));

nexttile
hold on
p = plot(log10(NFSD),log10(nfsd_vec),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
plot(log10(NFSD(1:12)),log10(roach_data_sep),'-o','MarkerFaceColor', 'k','Color', 'k','LineWidth',3)
plot(log10(NFSD(1:12)),log10(nfsd_vec2),'--s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);  
plot(log10(NFSD),log10(num_dafsd),'--s','LineWidth',3);  
%plot(log10(NFSD),log10(nfsd_vecai),':s','MarkerFaceColor', [0 0.4470 0.2],'LineWidth',3);
%plot(log10(NFSD),log10(nfsd_vecint),':o','MarkerFaceColor', [0 0.8 0.2],'LineWidth',10);
 legend({'$\int f^N(r,h) / g(h) dh$','Roach et al. (2018) Sep results', '$F^N(r_{12}) = \sum_{j = 12}^{16} F^N(r_j)$','$\int f^N(r,h) dh / \int g(h) dh$','aice','int'})
grid on
%title(sprintf("July"))
xticks(log10(xtick))
xticklabels(xticklab)
xlabel('Floe radius (m)')
yticks(log10(ytick))
yticklabels(yticklab)
%ylim([-9,13])
ylabel('Number distribution (km$^{-2}$)')
%xlim([0,3*10^3]);
hold off
f1.Position = [1500 800 600 500];
%exportgraphics(f1,'numdistlog.pdf','ContentType','vector')

f2 = figure(2);
map_plot(nfsd_cell,"aice",sector,grid,[0,600]);

%% Number FSD over time
sector = "SH";
grid = 'gx1'; 

case_name = 'wimoninit';
ssd_dir = '/Volumes/NoahDay5TB/cases/';
filedir = strcat(ssd_dir,case_name);
day = 01;
month_init = 7;
year = 2005;
date = sprintf('%d-0%d-0%d', year, month_init, day);
filename = strcat(filedir,"/history/iceh.",date,".nc");

data_2005jan = afsdn_to_nfsd(filename,sector);

ssd_dir = '/Volumes/NoahDay5TB/cases/';
filedir = strcat(ssd_dir,case_name);
day = 5;
month_init = 7;
year = 2005;
date = sprintf('%d-0%d-0%d', year, month_init, day);
filename = strcat(filedir,"/history/iceh.",date,".nc");

data_2005 = afsdn_to_nfsd(filename,sector);

day = 15;
month_init = 7;
year = 2005;
date = sprintf('%d-0%d-%d', year, month_init, day);
filename = strcat(filedir,"/history/iceh.",date,".nc");

data_2008 = afsdn_to_nfsd(filename,sector);

date = sprintf('%d-0%d-%d', 2005, 7, 30);
filename = strcat(filedir,"/history/iceh.",date,".nc");

data_2008dec = afsdn_to_nfsd(filename,sector);


roach_data_sep = [0.3*10^4, 3*10^1, 10^0, 0.8*10^(-1), 0.9*10^(-2), 1.2*10^(-3), 0.5*10^(-4), 0.7*10^(-5), 10^(-6), 1.3*10^(-7), 10^(-7), 10^(-4)];

f1 = figure(1);

t5 = tiledlayout(1,1);%(5,2);
t5.TileSpacing = 'compact';
cust_bounds =  max(NFSD);
xtick = 10.^(0:6);
xticklab = cellstr(num2str(round(log10(xtick(:))), '$10^{%d}$'));
ytick = 10.^(-20:2:13);
yticklab = cellstr(num2str(round(log10(ytick(:))), '$10^{%d}$'));

nfsd_vec2 = nfsd_vec(1:12-1);
nfsd_vec2(12) = sum(nfsd_vec(12:end));

nexttile
hold on
p = plot(log10(NFSD),log10(data_2005jan),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
plot(log10(NFSD(1:12)),log10(roach_data_sep),'-o','MarkerFaceColor', 'k','Color', 'k','LineWidth',3)
plot(log10(NFSD),log10(data_2005),'--s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);  
plot(log10(NFSD),log10(data_2008),'--s','LineWidth',3);  
plot(log10(NFSD),log10(data_2008dec),'--s','LineWidth',3);  
 legend({'01-07-2005','Roach et al. (2018) Sep results', '5-07-2005','15-07-2005','30-07-2005','int'})
grid on
xticks(log10(xtick))
xticklabels(xticklab)
xlabel('Floe radius (m)')
yticks(log10(ytick))
yticklabels(yticklab)
ylabel('Number distribution (km$^{-2}$)')
hold off
f1.Position = [1500 800 600 500];
exportgraphics(f1,'numdistlog.pdf','ContentType','vector')


%% 8. Time derivatives

% Get the NFSD from AFSDN
close all
clear xticks yticks f5 icemask f
aice_data = data_format_sector(filename,"aice",sector);
afsd = data_format_sector(filename,"afsd",sector);
icemask = aice_data > 0.01;
dafsd_latm_data(:,:,:) = data_format_sector(filename,"dafsd_latm",sector);
dafsd_latg_data(:,:,:) = data_format_sector(filename,"dafsd_latg",sector);
dafsd_newi_data(:,:,:) = data_format_sector(filename,"dafsd_newi",sector);
dafsd_weld_data(:,:,:) = data_format_sector(filename,"dafsd_weld",sector);
dafsd_wave_data(:,:,:) = data_format_sector(filename,"dafsd_wave",sector);

% dafsdrad_latm(:,:) = fsd_converter(filename,"dafsd_latm","fsdrad");
% dafsdrad_latg(:,:) = fsd_converter(filename,"dafsd_latg","fsdrad");
% dafsdrad_newi(:,:) = fsd_converter(filename,"dafsd_newi","fsdrad");
% dafsdrad_weld(:,:) = fsd_converter(filename,"dafsd_weld","fsdrad");
% dafsdrad_wave(:,:) = fsd_converter(filename,"dafsd_wave","fsdrad");

dafsdrad_latm(:,:) = fsd_converter(filename,"dafsd_latm","fsdrad");
dafsdrad_latg(:,:) = fsd_converter(filename,"dafsd_latg","fsdrad");
dafsdrad_newi(:,:) = fsd_converter(filename,"dafsd_newi","fsdrad");
dafsdrad_weld(:,:) = fsd_converter(filename,"dafsd_weld","fsdrad");
dafsdrad_wave(:,:) = fsd_converter(filename,"dafsd_wave","fsdrad");



raw_data = dafsd_newi_data(:,:,:);
for j = 1:384
    for i = 1:321
        work(i,j) = 0;
        for k = 1:Nf
             work(i,j) = work(i,j) + ...
                   raw_data(i,j,k);%*afsd(i,j,k);%*floe_binwidth(k)*(floe_rad_c(k)/aice(i,j));
        end
    end
end
dafsdrad_newi(:,:) = work;


dafsdrad_latm_mean = mean(dafsdrad_latm(icemask));
dafsdrad_latg_mean = mean(dafsdrad_latg(icemask));
dafsdrad_newi_mean = mean(dafsdrad_newi(icemask));
dafsdrad_weld_mean = mean(dafsdrad_weld(icemask));
dafsdrad_wave_mean = mean(dafsdrad_wave(icemask));
conFigure(11,3)
f1 = figure;
bar(1:5,[dafsdrad_latm_mean',dafsdrad_latg_mean',dafsdrad_newi_mean',dafsdrad_weld_mean',dafsdrad_wave_mean'])
xticks(1:5)
xticklabels({'latm','latg','newi','weld','wave'})
ylabel('$dr_a/dt$ (m/day)')
%legend({'latm','latg','newi','weld','wave'},'Location','north')
%f1.Position = [1200 100 600 500];
%exportgraphics(f1,'numdistlog.pdf','ContentType','vector')

max_lim = 10;
f = figure;
t1 = tiledlayout(1,5);
nexttile
max_lim = max(max(dafsdrad_latm));
if max_lim == 0 
    max_lim = eps;
end
map_plot(dafsdrad_latm,"dafsd_latm",sector,grid,[-max_lim,max_lim]);
cmocean('balance')
title("Lat melt")
nexttile
max_lim = max(max(dafsdrad_latg));
map_plot(dafsdrad_latg,"dafsd_latg",sector,grid,[-max_lim,max_lim]);
cmocean('balance')
title("Lat growth")
nexttile
max_lim = max(max(dafsdrad_newi));
map_plot(dafsdrad_newi,"dafsd_newi",sector,grid,[-max_lim,max_lim]);
cmocean('balance')
title("New ice")
nexttile
max_lim = max(max(dafsdrad_weld));
map_plot(dafsdrad_weld,"dafsd_weld",sector,grid,[-max_lim,max_lim]);
cmocean('balance')
title("Weld")
nexttile
max_lim = max(max(abs(dafsdrad_wave)))/2;
if max_lim == 0
    max_lim = eps;
end
map_plot(dafsdrad_wave,"dafsd_wave",sector,grid,[-max_lim,max_lim]);
cmocean('balance')
title("Wave")
nexttile
%%
aice = data_format_sector(filename,"aice",sector);
afsd2 = afsd + dafsd_latm_data + dafsd_latg_data + dafsd_newi_data + dafsd_weld_data + dafsd_wave_data;
raw_data = afsd2(:,:,:);
for j = 1:384
    for i = 1:321
        work(i,j) = 0;
        for k = 1:Nf
             work(i,j) = work(i,j) + ...
                   raw_data(i,j,k)*floe_binwidth(k)*(floe_rad_c(k)/aice(i,j));
        end
    end
end
fsdrad_change(:,:) = work;
f = figure;
map_plot(fsdrad_change,"fsdrad",sector,grid,[0,10000]);


%%
num_dafsd_latm = numberfsdconverter(filename,dafsd_latm_data);
num_dafsd_latg = numberfsdconverter(filename,dafsd_latg_data);
num_dafsd_newi = numberfsdconverter(filename,dafsd_newi_data);
num_dafsd_weld = numberfsdconverter(filename,dafsd_weld_data);
num_dafsd_wave = numberfsdconverter(filename,dafsd_wave_data);

plot(1:16,[num_dafsd_latm',num_dafsd_latg',num_dafsd_newi',num_dafsd_weld',num_dafsd_wave'])

%% 9. Perimeter density comparison with Bateson et al. THIS DOESNT WORK, CICE DOESNT OUTPUT THE PERIMETER FSD
% Get the NFSD from AFSDN
%close all
clear xticks yticks f5 icemask nfsd trcrn
aice_data = data_format_sector(filename,"aice",sector);
icemask0 = aice_data > -eps;
icemask = aice_data > 0.01;
icemask90 = aice_data > 0.9;
fsdperim_data(:,:,:,:) = data_format_sector(filename,"fsdperim",sector);
afsd_data(:,:,:) = data_format_sector(filename,"afsd",sector);
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
tarea = data_format(filename,"tarea");

[nx,ny,nf,nc] = size(afsdn_data);
%icemask = icemask0;


% Make a mask for the Southern Hemisphere
% Coordinates
ocean_mask = data_format(filename,'tmask');
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
sector_mask = false(nx,ny);
for i = 0:lat_out(1)-lat_out(2) % Cycle through latitudes
    sector_mask(lon_out(1):lon_out(3), i+1) = true;
end
sector_mask = logical(sector_mask.*ocean_mask);




% for k = 1:nf
%     temp = fsdperim_data(:,:,k);
%     nfsd_vec(k) = mean(temp(sector_mask)); % Take the cells over the sector
%     temp = nfsd2(:,:,k);
%     nfsd_vecai(k) = mean(temp(sector_mask)); % Take the cells over the sector
%     temp = nfsd3(:,:,k);
%     nfsd_vecint(k) = mean(temp(sector_mask)); % Take the cells over the sector
% end
%
[N,edges] = histcounts(fsdperim_data);
%num_dafsd = numberfsdconverter(filename,afsd_data);

%roach_data_sep = [0.3*10^4, 3*10^1, 10^0, 0.8*10^(-1), 0.9*10^(-2), 1.2*10^(-3), 0.5*10^(-4), 0.7*10^(-5), 10^(-6), 1.3*10^(-7), 10^(-7), 10^(-4)];

f1 = figure(1);

t5 = tiledlayout(1,1);%(5,2);
t5.TileSpacing = 'compact';
cust_bounds =  max(NFSD);
xtick = 10.^(0:6);
xticklab = cellstr(num2str(round(log10(xtick(:))), '$10^{%d}$'));
ytick = 10.^(-9:2:13);
yticklab = cellstr(num2str(round(log10(ytick(:))), '$10^{%d}$'));

conFigure(11)
p = plot(log10(edges(2:end)),log10(N.*10^3),'-s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);
%p = bar(log10(edges(2:end)),log10(N.*10^3))
%plot(log10(NFSD(1:12)),log10(roach_data_sep),'-o','MarkerFaceColor', 'k','Color', 'k','LineWidth',3)
%plot(log10(NFSD(1:12)),log10(nfsd_vec2),'--s','MarkerFaceColor', [0 0.4470 0.7410],'LineWidth',3);  
%plot(log10(NFSD),log10(num_dafsd),'--s','LineWidth',3);  
%plot(log10(NFSD),log10(nfsd_vecai),':s','MarkerFaceColor', [0 0.4470 0.2],'LineWidth',3);
%plot(log10(NFSD),log10(nfsd_vecint),':o','MarkerFaceColor', [0 0.8 0.2],'LineWidth',10);
 legend({'September','aice','int'})
grid on
%title(sprintf("July"))
xticks(log10(xtick))
xticklabels(xticklab)
xlabel('Floe radius (m)')
yticks(log10(ytick))
yticklabels(yticklab)
%ylim([-9,13])
ylabel('Perimenter density (m$^{-1}$)')
%xlim([0,3*10^3]);
hold off
%f1.Position = [1500 800 600 500];
%exportgraphics(f1,'numdistlog.pdf','ContentType','vector')













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

function y = stat_vec(data)
% ["Mean","Median","Mode","Standard deviation","Maximum"];
    y = [mean(data); median(data); mode(data); std(data); max(data)];
end

% -------------------------------------------------------------------------

function output = numberfsdconverter(filename,data)
NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
[lat, lon, row] = grid_read("gx1");
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;
aice_data = data_format_sector(filename,"aice","SH");
icemask = aice_data > 0.01;
[nx,ny,nf] = size(data);
% Make a mask for the Southern Hemisphere
% Coordinates
sector = "SH";
ocean_mask = data_format(filename,'tmask');
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
sector_mask = false(nx,ny);
for i = 0:lat_out(1)-lat_out(2) % Cycle through latitudes
    sector_mask(lon_out(1):lon_out(3), i+1) = true;
end
sector_mask = logical(sector_mask.*ocean_mask);
% Get tracer array
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            if aice_data(i,j) < eps
                trcrn(i,j,k) = 0;
            else
                trcrn(i,j,k) = data(i,j,k)*floe_binwidth(k)*(10^(-3))/aice_data(i,j); % something x Length: (something m)
            end
        end
    end
end

% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*(floe_rad_c*10^(-3)).^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
                nfsd(i,j,k) = trcrn(i,j,k)/floe_area_c(k); % Dimensionless: something / Length
        end
    end
end

% Intergrate nfstd to get nfsd

for k = 1:nf
    temp = nfsd(:,:,k);
    nfsd_vec(k) = mean(temp(sector_mask));
end
output = nfsd_vec;
end

% -------------------------------------------------------------------------
function [out_data] = afsdn_to_nfsd(filename,sector)
NFSD = ncread(filename,"NFSD");
NCAT = ncread(filename,"NCAT");
[lat, lon, row] = grid_read("gx1");
Nf = 16;
Nc = 5;
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;
aice_data = data_format_sector(filename,"aice",sector);
icemask = aice_data > 0.01;9;
afsdn_data(:,:,:,:) = data_format_sector(filename,"afsdn",sector);
afsd_data(:,:,:) = data_format_sector(filename,"afsd",sector);
aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
tarea = data_format(filename,"tarea");

[nx,ny,nf,nc] = size(afsdn_data);


% Make a mask for the Southern Hemisphere
% Coordinates
ocean_mask = data_format(filename,'tmask');
coords = sector_coords(sector); % (NW;NE;SW;SW) (lat,lon) %(NW;SW,NE,SE)
for i = 1:4
    [lat_out(i),lon_out(i)] = lat_lon_finder(coords(i,1),coords(i,2),lat,lon);
end
sector_mask = false(nx,ny);
for i = 0:lat_out(1)-lat_out(2) % Cycle through latitudes
    sector_mask(lon_out(1):lon_out(3), i+1) = true;
end
sector_mask = logical(sector_mask.*ocean_mask);

% Get tracer array
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            for n = 1:nc
                if aicen_data(i,j,n) < eps
                    trcrn(i,j,k,n) = 0;
                else
                    trcrn(i,j,k,n) = afsdn_data(i,j,k,n)*floe_binwidth(k)*(10^(-3))/aicen_data(i,j,n); % something x Length: (something m) now oonverted to km
                end
            end
        end
    end
end

% Calculate nfstd
alpha = 0.66; % Dimensionless
floe_area_c = 4*alpha*(floe_rad_c*(10^(-3))).^2; % Area: (m^2)
for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            for n = 1:nc
                % Get the number of floes per cell
                nfstd(i,j,k,n) = trcrn(i,j,k,n)/(floe_area_c(k)); % Dimensionless: something / Length
            end
        end
    end
end

% Intergrate nfstd to get nfsd

for i = 1:nx
    for j = 1:ny
        for k = 1:nf
            nfsd(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
        end
         nfsd_cell(i,j) = sum(nfsd(i,j,:));
    end
end

for k = 1:nf
    temp = nfsd(:,:,k);
    nfsd_vec(k) = mean(temp(sector_mask)); % Take the cells over the sector
end
 out_data = nfsd_vec;
end