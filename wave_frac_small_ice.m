% Plot the amounts of ice being put into the smallest FSD category from the
% wave breakup.

clear all
clc
close all
addpath functions

historydir = '/Users/noahday/GitHub/CICE-plotting-tools/cases/monthwim/history/';%'/Volumes/NoahDay5TB/cases/wimoninit/history/';

a = dir([historydir '/*.nc']);
n_files = numel(a);

for i = 1:n_files
   filenames(i,:) = strcat(historydir,a(i).name);
   dirdates(i,:) = a(i).name(6:end-3);
end