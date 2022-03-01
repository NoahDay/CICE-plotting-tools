function label = colorlabel(variable)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if variable == "aice"
    label = "Concentration";
elseif variable == "vatm"
    label = "Wind speed (m/s)";
elseif variable == "wave_sig_ht"
    label = "SWH (m)";
elseif variable == "fsdrad"
    label = "Floe size radius with SIC > 0.01 (m)";
elseif variable == "Tair" || variable == "Tsfc" || variable == "sst"
    label = "Degrees C";
elseif variable == "vvel"
    label = "Ice drift speed (m/s)";
else
    label = "undefined";
end
end