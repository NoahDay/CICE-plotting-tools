function max_wind = storm_finder(filename, grid, sector)
%STORM_FINDER Determines if a storm is present (by wind velocity)

var1 = "uatm"; % atm velocity (x)
var2 = "vatm"; % atm velocity (y)
dim = 2;
atm_data_x = data_format_sector(filename,var1,sector);
atm_data_y = data_format_sector(filename,var2,sector);
magnitude_wind = sqrt(atm_data_x.^2+atm_data_y.^2);
max_wind = max(magnitude_wind(1:numel(magnitude_wind)));
end