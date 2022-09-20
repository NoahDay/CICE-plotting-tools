% Plot the stress vector diagram for a cell
% 1. Preamble
% 2. Read in data
%   2a. Ice drift
%   2b. Atm stress
%   2c. Ocn stress
%   2d. Cor stress
%   2e. Sea surface slope tlt stress
%   2f. Internal stress
% 3. Plot the vector diagram

%% 1. Preamble
close all
clear all
addpath functions

user = 'noahday'; %a1724548, noahday, Noah
case_name = 'ocnforcing';
grid = 'gx1'; 
day = 3;
month = 7;
year = 2009;
sector = "SA";
if day < 9
    date = sprintf('%d-0%d-0%d', year, month, day);
else
    date = sprintf('%d-0%d-%d', year, month, day);
end

dim = 2;
[lat,lon,row] = grid_read(grid);


ssd = 1;
if ssd == 1
    filename = strcat('/Volumes/NoahDay5TB/cases/ocnforcing/history/iceh.',date,".nc");
else
    filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
end

aice_data = data_format_sector(filename,"aice",sector);
lat_vec = reshape(lat,1,[]);
lon_vec = reshape(lon,1,[]);
SIC = 0.15;
[lat_ice_edge, lon_ice_edge, edge] = find_ice_edge(aice_data,SIC,sector,lat,lon);

lon_pos = 30;
%% 2. Read in data

% 2a. Ice drift
Data.uvel = data_format(filename,"uvel",row,lat,lon,dim);
Data.vvel = data_format(filename,"vvel",row,lat,lon,dim);

%   2b. Atm stress
Data.strairx = data_format(filename,"strairx",row,lat,lon,dim);
Data.strairy = data_format(filename,"strairy",row,lat,lon,dim);

%   2c. Ocn stress
Data.strocnx = data_format(filename,"strocnx",row,lat,lon,dim);
Data.strocny = data_format(filename,"strocny",row,lat,lon,dim);

%   2d. Cor stress
Data.strcorx = data_format(filename,"strcorx",row,lat,lon,dim);
Data.strcory = data_format(filename,"strcory",row,lat,lon,dim);

%   2e. Sea surface slope tlt stress
Data.strtltx = data_format(filename,"strtltx",row,lat,lon,dim);
Data.strtlty = data_format(filename,"strtlty",row,lat,lon,dim);

%   2f. Internal stress
Data.strintx = data_format(filename,"strintx",row,lat,lon,dim);
Data.strinty = data_format(filename,"strinty",row,lat,lon,dim);

%% 3. Plot the vector diagram
% lon(30,:) ~ 32.5

% Data at (-58 S, 32.5 E)
% cell_uvel = [0, data_uvel(30,edge(30))];
% cell_vvel = [0,data_vvel(30,edge(30))];
% 
% %   2b. Atm stress
% cell_strairx = data_strairx(30,edge(30));
% cell_strairy = data_strairy(30,edge(30));
% 
% %   2c. Ocn stress
% cell_strocnx = data_strocnx(30,edge(30));
% cell_strocny = data_strocny(30,edge(30));
% 
% %   2d. Cor stress
% cell_strcorx = data_strcorx(30,edge(30));
% cell_strcory = data_strcory(30,edge(30));
% 
% %   2e. Sea surface slope tlt stress
% cell_strtltx = data_strtltx(30,edge(30));
% cell_strtlty = data_strtlty(30,edge(30));
% 
% %   2f. Internal stress
% cell_strintx = data_strintx(30,edge(30));
% cell_strinty = data_strinty(30,edge(30));
edge_vec = edge(30)-1:-1:24;
n = length(edge_vec);
t = tiledlayout(2,ceil(n/2));
line_width = 2;
for i = edge_vec
    lat_pos = i;
    vectors = get_stresses(Data, lon_pos, lat_pos);
    nexttile
    stress_plot(vectors)
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    t = title(sprintf("(%g S, %g E)", lat(lon_pos,lat_pos), lon(lon_pos,lat_pos)));
    t.FontSize = 10;
end

%% Functions
function vectors = get_stresses(data, lon_pos, lat_pos)
    % Get the stress data
   % cell_uvel = [0, data_uvel(30,edge(30))];
   % cell_vvel = [0,data_vvel(30,edge(30))];
    
    %   2b. Atm stress
    cell_strairx = data.strairx(lon_pos,lat_pos);
    cell_strairy = data.strairy(lon_pos,lat_pos);
    
    %   2c. Ocn stress
    cell_strocnx = data.strocnx(lon_pos,lat_pos);
    cell_strocny = data.strocny(lon_pos,lat_pos);
    
    %   2d. Cor stress
    cell_strcorx = data.strcorx(lon_pos,lat_pos);
    cell_strcory = data.strcory(lon_pos,lat_pos);
    
    %   2e. Sea surface slope tlt stress
    cell_strtltx = data.strtltx(lon_pos,lat_pos);
    cell_strtlty = data.strtlty(lon_pos,lat_pos);
    
    %   2f. Internal stress
    cell_strintx = data.strintx(lon_pos,lat_pos);
    cell_strinty = data.strinty(lon_pos,lat_pos);

    vectors = [cell_strairx, cell_strocnx, cell_strcorx, cell_strtltx, cell_strintx;
        cell_strairy, cell_strocny, cell_strcory, cell_strtlty, cell_strinty];
end

function stress_plot(vectors)
    line_width = 2;
    quiver(0,0,vectors(1,1),vectors(2,1),'linewidth',line_width)
    hold on
    quiver(0,0,vectors(1,2),vectors(2,2),'linewidth',line_width)
    quiver(0,0,vectors(1,3),vectors(2,3),'linewidth',line_width)
    quiver(0,0,vectors(1,4),vectors(2,4),'linewidth',line_width)
    quiver(0,0,vectors(1,5),vectors(2,5),'linewidth',line_width)
        set(gcf,'Position',[1000 1000 500 400])
    hold off
    legend({"Air stress", "Ocean stress", "Coriolis stress", "Sea surface slope stress", "Internal stress"})
end