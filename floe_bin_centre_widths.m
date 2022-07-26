%% Set up
clc; clear all; close all
filename = '/Users/noahday/GitHub/CICE-plotting-tools/cases/pancake_tracer/history/iceh.2005-01.nc';
addpath functions/
NCAT = ncread(filename,"NCAT"); NFSD = ncread(filename,"NFSD");
floe_rad_c = NCAT;

Nf = numel(NFSD); % 'category floe size (center)'
Nc = numel(NCAT); % 'category maximum thickness'

% Set upc
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2; % = NFSD


% NFSD == floe_rad_c
%figure
%bar([NFSD'-floe_rad_c])


min_floe = 6.65000000e-02;
fsd_lims(1) = min_floe;
for i = 1:numel(NFSD)
    fsd_lims(i+1) = 2*NFSD(i)-fsd_lims(i);
end