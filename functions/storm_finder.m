function [max_wind, mean_wind, min_wind] = storm_finder(filename, sector)
%STORM_FINDER Determines if a storm is present within a sector (by wind velocity)

    var1 = "uatm"; % atm velocity (x)
    var2 = "vatm"; % atm velocity (y)
    dim = 2;
    atm_data_x = data_format_sector(filename,var1,sector);
    atm_data_y = data_format_sector(filename,var2,sector);
    magnitude_wind = sqrt(atm_data_x.^2+atm_data_y.^2);
    max_wind = max(magnitude_wind(1:numel(magnitude_wind)));
    idx = isnan(magnitude_wind);
    mean_wind = mean(magnitude_wind(~idx));
    min_wind = min(magnitude_wind(1:numel(magnitude_wind)));
end