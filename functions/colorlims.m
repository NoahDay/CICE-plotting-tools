function limit = colorlims(variable)
%COLORLIMS given a variable output the appropriate colorbar limits.
% Author: Noah Day 30 Sep 21

 if variable == "wave_sig_ht"
    limit=[0,6];
 elseif variable == "wave_sig_ht_d"
    limit=[0,6];
 elseif variable == "fsdrad"
    limit=[0,1000];
 elseif variable ==  "fsdrad_d"
    limit=[0,2000];
 elseif variable  == "aice"
    limit=[0,1];
 elseif variable  == "aice_d"
    limit=[0,1];
 elseif variable  == "peak_period"
    limit=[0,20];
elseif variable  == "peak_period_d"
    limit=[0,20];
 elseif variable == "mean_wave_dir"
     limit=[0,2*pi];
 elseif variable == "mean_wave_dir_d"
     limit=[0,2*pi];
 elseif variable == "hi"
     limit = [0,2];
 elseif variable == "frazil"
     limit = [0,5];
 elseif variable == "sigP"
     limit = [0,5000];
 elseif variable == "Tair"
     limit = [-30,30];
 elseif variable == "sst"
     limit = [-10,15];
 elseif variable == "Tsfc"
     limit = [-20,4];
 elseif variable == "dafsd_wave"
     limit = [-0.05,0.05];
 elseif variable == "iage"
     limit = [0,2];
 elseif variable == "vatm"
     limit = [0,35];
 else
     limit = [0,1];
 end
end

