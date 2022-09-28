filename = '/Users/noahday/Downloads/IBTrACS.NI.v04r00.nc';
ncdisp(filename)

lat = ncread(filename,'lat'); % Degrees N
lon = ncread(filename,'lon'); % Degrees E