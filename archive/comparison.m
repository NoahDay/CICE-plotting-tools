clear all
close all

% Comparing CICE results with and without waves
 plot_title1 = 'benchmark';
 plot_title2 = 'wim on';
 case1dir = '1yearnowaves';
 case2dir = '8month';
 date = '2005-06-19';
 datapoints = 1;
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'noahday'; %a1724548, noahday, Noah
 variable = 'wave_sig_ht'; % wave_sig_ht, peak_period, fsdrad, 
 map_type = 'eqaazim'; %cassini
 rad = 110;
 
 % sea ice area: 'aicen' or 'aice'
 % horizonal velocity: 'uvel'
 % vertical velocity: 'vvel'
 
filedir1 = strcat('cases/',case1dir,'/history/iceh.',date,'.nc');
filedir2 = strcat('cases/',case2dir,'/history/iceh.',date,'.nc');

f1 = figure;
comparison_maker(filedir1, filedir2, plot_title1, plot_title2, variable, grid, map_type, user);

f2 = figure;
difference_maker(filedir1, filedir2, plot_title1, variable, grid, map_type, user)

plot_title_pan = sprintf('Floes with radius smaller than %d m',rad); 

if variable == "fsdrad"
    f3 = figure;
    pancake_maker(filedir1, plot_title_pan, variable, grid, map_type, user, rad)
end