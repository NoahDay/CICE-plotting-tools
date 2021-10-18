function limit = colorlims(variable)
%COLORLIMS given a variable output the appropriate colorbar limits.
% Author: Noah Day 30 Sep 21

 if variable == "wave_sig_ht"
    limit=[0,10];
 elseif variable == "wave_sig_ht_d"
    limit=[0,10];
 elseif variable == "fsdrad"
    limit=[0,400];
 elseif variable ==  "fsdrad_d"
    limit=[0,800];
 elseif variable  == "aice"
    limit=[0,1];
 elseif variable  == "aice_d"
    limit=[0,1];
 elseif variable  == "peak_period"
    limit=[0,20];
 elseif variable == "mean_wave_dir"
     limit=[0,2*pi];
 else
     limit = [0,1];
 end
end

