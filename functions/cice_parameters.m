function [floe_binwidth, floe_rad_l, floe_rad_h] = cice_parameters(NFSD)
%CICE_PARAMETERS - Given a FSD bins centres vector will return the binwidth
%
% Syntax:  [floe_binwidth, floe_rad_l, floe_rad_h] = cice_parameters(NFSD)
%
% Inputs:
%    NFSD - floe bin centres, vector (1xn) [m]
%
% Outputs:
%    floe_binwidth - floe bin widths, vector (1xn) [m]
%    floe_rad_l    - floe bin lower lims, vector (1xn) [m]
%    floe_rad_h    - floe bin upper lims, vector (1xn) [m]
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Noah Day
% Email: noah.day@adelaide.edu.au
% July 2022; Last revision: 26-July-2022

%------------- BEGIN CODE --------------

min_floe = 6.65000000e-02;
fsd_lims(1) = min_floe;
for i = 1:numel(NFSD)
    fsd_lims(i+1) = 2*NFSD(i)-fsd_lims(i);
end

floe_rad_l = fsd_lims(1:end-1);
floe_rad_h = fsd_lims(2:end);
floe_binwidth = floe_rad_h - floe_rad_l;
%------------- END OF CODE --------------

end