%% Determining the MIZ
% The MIZ is calculated from a combination of FSD radius, significant wave
% height, and the change in FSD due to waves.
clear all
close all
addpath functions
addpath packages/bedmap2_toolbox_v4
% Switches
sic_miz_switch = 0;
swh_miz_switch = 0;
fsd_miz_switch = 0;
miz_width_switch = 1;
iceedge = 0;

plotting = 0;
plotting_mean = 1;
printing = 0;
sector = [-65, 360-132.2];
ssd = 0; % ssd or local data
SIC = [0.15, 0.8];
SWH = 0.0001;
fsd_max = 10; % Maximum size of pancake ice
ymax = 1200;
%% Preamble
user = 'noahday'; %a1724548, noahday, Noah
case_name = 'fixedwaves';
if ssd == 1
    ssd_dir = '/Volumes/Noah_SSD/run_data';
    filedir = strcat(ssd_dir,case_name);
else
    filedir = strcat('cases/',case_name);
end
grid = 'gx1'; 
time_period = 'd'; %'1','d','m','y'
season = "winter";

day = 1;
month_init = 1;
year = 2006;
date = sprintf('%d-0%d-0%d', year, month_init, day);
months = [2,5,9,12];
dim = 2;
% Load the grid, ice shelf, and MIZ statistics.
% Grid
[lat,lon,row] = grid_read(grid);

% % Load ice shelf data
% shelf = bedmap2_data('icemask');
% shelf = 0.0*(shelf == 1); 
% shelf(shelf==0) = NaN; % 1 iceshelf
% [latshelf,lonshelf] = bedmap2_data('latlon');
figcount = 0;


%% 1. SIC defintion
% The SIC MIZ is defined between 0.15 and 0.8 SIC
if sic_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 2;
    variable = "aice";
    sic_data = data_format_sector(filename,variable,"SA");
    idx = sic_data < SIC(1); 
    sic_miz = sic_data;
    sic_miz(idx) = 0.0;
    idx = sic_data > SIC(2); 
    sic_miz(idx) = 0.0;
    figcount = figcount + 1;
    figure(figcount)
    map_plot(sic_data,variable,"SA",grid);
end

%% 2. SWH definition
if swh_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 2;
    variable = "wave_sig_ht";
    swh_data = data_format_sector(filename,variable,sector);
    idx = swh_data < eps; 
    swh_miz = swh_data;
    swh_miz(idx) = 0.0;
    figcount = figcount + 1;
    figure(figcount)
    map_plot(swh_miz,variable,sector,grid);
end

%% 3. FSD definition
figcount = 0;
if fsd_miz_switch == 1
    date = sprintf('%d-0%d-0%d', year, month_init, day);
    filename = strcat(filedir,"/history/iceh.",date,".nc");
    [len, wid] = size(lat);
    dim = 3;
    variable = "afsd";
    NFSD = ncread(filename,"NFSD");
    fsd_data = data_format_sector(filename,variable,sector,dim);
%     idx = fsd_data > 10; 
%     fsd_miz = fsd_data;
%     fsd_miz(idx) = 0.0;
    for i = 1:len
        for j = 1:wid
            fsd_pdf = zeros(1,length(NFSD));
            for k = 1:length(NFSD)
                fsd_pdf(k) = fsd_data(i,j,k);
            end
            fsd_miz(i,j) = fsd_pdf*(NFSD.^2*pi); % Fraction x Area
            %if ~isnan(fsd_pdf*(NFSD.^2*pi))
            %    temp = fsd_pdf.*(NFSD.^2*pi);
            %end
        end
    end
    figcount = figcount + 1;
    figure(figcount)
    map_plot(fsd_miz,variable,sector,grid,[0,1000]);

    % Normalize FSD per cell
    variable = "afsd";
    NFSD = ncread(filename,"NFSD");
    NCAT = ncread(filename,"NCAT");
    dim = 3;
    fstd_data = data_format_sector(filename,variable,sector,dim); % (x,y,fsd,itd)
    for i = 1:len
        for j = 1:wid
            fsd_pdf = zeros(1,length(NFSD));
            for k = 1:length(NFSD)
%                 % Normalize ITD
%                 for l = 1:length(NCAT)
%                     itd_pdf(l) = fstd_data(i,j,k,l)/sum(fstd_data(i,j,k,:));
%                 end
%                 % Integrate the ITD
%                 Int_itd = itd_pdf*NCAT;
%                 fsd_data = fstd_data(:,:,:,1);
                % Intergrate w.r.t. FSD
                fsd_pdf(k) = fsd_data(i,j,k)/sum(fsd_data(i,j,:));
            end
            fsd_miz(i,j) = fsd_pdf*NFSD; % Weighted Average
        end
    end
    figcount = figcount + 1;
    figure(figcount)
    map_plot(fsd_miz,variable,sector,grid,[0,1000]);
    
    
    % FSDrad 
    variable = "fsdrad";
    fsdrad_data = data_format_sector(filename,variable,sector,2);
    figcount = figcount + 1;
    figure(figcount)
    map_plot(fsdrad_data,variable,sector,grid,[0,1000]);

    % AICE 
    variable = "aice";
    aice_data = data_format_sector(filename,variable,sector,2);
    figcount = figcount + 1;
    figure(figcount)
    map_plot(aice_data,variable,sector,grid,[0,1]);

    figcount = figcount + 1;
    figure(figcount)
    map_plot(fsdrad_data-fsd_miz,variable,sector,grid,[0,1000]);
end
%% 4. MIZ widths
if miz_width_switch == 1
    Data = struct('Month',{},'SIC',{},'SWH',{},'FSD',{});
    % Sector analysis:
    if isstring(sector) 
        coords = sector_coords(sector);
        for k = 1:length(months)
        month_init = months(k);
        if month_init < 10
            date = sprintf('%d-0%d-0%d', year, month_init, day);
        else
            date = sprintf('%d-%d-0%d', year, month_init, day);
        end
        if month_init == 1 || month_init == 3 || month_init == 5 || month_init == 7 || month_init == 8 || month_init == 10 ||  month_init == 12
            datapoints = 31;
        elseif month_init == 2
            datapoints = 28;
        else
            datapoints = 30;
        end
        distSIC = zeros(1,datapoints);
        distSWH = zeros(1,datapoints); 
        distFSD = zeros(1,datapoints);
        for j = 1:datapoints
            filename = strcat(filedir,"/history/iceh.",date,".nc");
           
            % Find the longitudes of the edges of the sector
            [~,lon_out1] = lat_lon_finder(coords(1,1),coords(1,2),lat,lon); 
            [~,lon_out2] = lat_lon_finder(coords(3,1),coords(3,2),lat,lon); 
            
            
            % 4. a) SIC
            variable = "aice";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
            transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
            
    
            for l = 1:lon_out2-lon_out1
                % Find the cells that satisfy the requirements
                idx1 = transect_data(l,:) > SIC(1);
                idx2 = transect_data(l,:) < SIC(2);
                idx_transect = logical(idx1.*idx2);
                
                points = lat(lon_out1+l-1,idx_transect);
                southern_hemi_points = points(points < 0);
                [~,wid] = size(southern_hemi_points);
                if wid == 0 % No MIZ
                    distSIC(l,j) = 0;
                else
                    distSIC(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                end
                if printing == 1
                    fprintf(strcat("The width of the SIC MIZ in the ", sector," sector is: %g km\n"), distSIC)
                end
            end
            
            % 4.b SWH
            variable = "wave_sig_ht";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
            transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
            
    
            for l = 1:lon_out2-lon_out1
                % Find the cells that satisfy the requirements
               idx_transect = transect_data(l,:) > SWH;
                
                points = lat(lon_out1+l-1,idx_transect);
                southern_hemi_points = points(points < 0);
                [~,wid] = size(southern_hemi_points);
                if wid == 0 % No MIZ
                    distSWH(l,j) = 0;
                else
                    distSWH(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                end
                if printing == 1
                    fprintf(strcat("The width of the SWH MIZ in the ", sector," sector is: %g km\n"), distSWH)
                end
            end
    
    
    
            % 4.c FSD
            variable = "fsdrad";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        
            end
            transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
            
    
            for l = 1:lon_out2-lon_out1
                % Find the cells that satisfy the requirements
                idx1 = transect_data(l,:) > eps;
                idx2 = transect_data(l,:) < fsd_max;
                idx_transect = logical(idx1.*idx2);
                
                points = lat(lon_out1+l-1,idx_transect);
                southern_hemi_points = points(points < 0);
                [~,wid] = size(southern_hemi_points);
                if wid == 0 % No MIZ
                    distFSD(l,j) = 0;
                else
                    distFSD(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                end
                if printing == 1
                    fprintf(strcat("The width of the FSD MIZ in the ", sector," sector is: %g km\n"), distFSD)
                end
            end
            date = update_date(date);
        end
        Data(k).distSIC = distSIC;
        Data(k).distSWH = distSWH;
        Data(k).distFSD = distFSD;
        Data(k).Month = monthName(months(k));
        disp('complete')
    end
    else % is a transect
        for k = 1:length(months)
        month_init = months(k);
        if month_init < 10
            date = sprintf('%d-0%d-0%d', year, month_init, day);
        else
            date = sprintf('%d-%d-0%d', year, month_init, day);
        end
        if month_init == 1 || month_init == 3 || month_init == 5 || month_init == 7 || month_init == 8 || month_init == 10 ||  month_init == 12
            datapoints = 31;
        elseif month_init == 2
            datapoints = 28;
        else
            datapoints = 30;
        end
        datapoints = 22;
        distSIC = zeros(1,datapoints);
        distSWH = zeros(1,datapoints); 
        distFSD = zeros(1,datapoints);
        for j = 1:datapoints
            filename = strcat(filedir,"/history/iceh.",date,".nc");
            % 4. a) SIC
            [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "aice";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
            transect_data = output_data(lon_out,:); % Take the southern hemisphere
            
            
            idx1 = transect_data > SIC(1);
            idx2 = transect_data < SIC(2);
            idx_transect = logical(idx1.*idx2);
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distSIC(j) = 0;
            else
                distSIC(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the SIC MIZ along %g E is %g km\n', sector(2), distSIC)
            end
            
            % 4.b SWH
            [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "wave_sig_ht";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
        
            transect_data = output_data(lon_out,:); % Take the southern hemisphere
            
            idx_transect = transect_data > SWH;
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distSWH(j) = 0;
            else
                distSWH(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the SWH MIZ along %g E is %g km\n', sector(2), distSWH)
            end
            % 4.c FSD
            
            [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "fsdrad";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        
            end
            transect_data = output_data(lon_out,:); % Take the southern hemisphere
            idx1 = transect_data < fsd_max;
            idx2 = transect_data > eps;
            idx_transect = logical(idx1.*idx2);
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distFSD(j) = 0;
            else
                distFSD(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the FSD MIZ along %g E is %g km\n', sector(2), distFSD)
            end
            date = update_date(date);
        end
        Data(k).SIC = distSIC;
        Data(k).SWH = distSWH;
        Data(k).FSD = distFSD;
        Data(k).Month = monthName(months(k));
    end
    end
end

%% 5. MIZ widths with SIC
if miz_width_switch == 1
    Data_sic = struct('Month',{},'SIC',{},'SWH',{},'FSD',{});
    % Sector data
    if isstring(sector) 
        coords = sector_coords(sector);
        for k = 1:length(months)
            month_init = months(k);
            if month_init < 10
                date = sprintf('%d-0%d-0%d', year, month_init, day);
            else
                date = sprintf('%d-%d-0%d', year, month_init, day);
            end
            if month_init == 1 || month_init == 3 || month_init == 5 || month_init == 7 || month_init == 8 || month_init == 10 ||  month_init == 12
                datapoints = 31;
            elseif month_init == 2
                datapoints = 28;
            else
                datapoints = 30;
            end
    %         datapoints = 22;
            distSIC = zeros(1,datapoints);
            distSWH = zeros(1,datapoints); 
            distFSD = zeros(1,datapoints);
            for j = 1:datapoints
                filename = strcat(filedir,"/history/iceh.",date,".nc");
               
                % Find the longitudes of the edges of the sector
                [~,lon_out1] = lat_lon_finder(coords(1,1),coords(1,2),lat,lon); 
                [~,lon_out2] = lat_lon_finder(coords(3,1),coords(3,2),lat,lon); 
                
                
                % 4. a) SIC
                variable = "aice";
                data = data_format(filename,variable,row,lat,lon,dim);
                if plotting == 1
                    figcount = figcount + 1;
                    figure(figcount)
                    [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                        title('FSDrad along transect')
                else
                    [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                end
                transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
                
        
                for l = 1:lon_out2-lon_out1
                    % Find the cells that satisfy the requirements
                    idx1 = transect_data(l,:) > SIC(1);
                    idx2 = transect_data(l,:) < SIC(2);
                    idx_transect = logical(idx1.*idx2);
                    
                    points = lat(lon_out1+l-1,idx_transect);
                    southern_hemi_points = points(points < 0);
                    [~,wid] = size(southern_hemi_points);
                    if wid == 0 % No MIZ
                        distSIC(l,j) = 0;
                    else
                        distSIC(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                    end
                    if printing == 1
                        fprintf(strcat("The width of the SIC MIZ in the ", sector," sector is: %g km\n"), distSIC)
                    end
                end
                
                % 4.b SWH
                variable = "wave_sig_ht";
                data = data_format(filename,variable,row,lat,lon,dim);
                if plotting == 1
                    figcount = figcount + 1;
                    figure(figcount)
                    [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                else
                    [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                end
                transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
                
        
                for l = 1:lon_out2-lon_out1
                    % Find the cells that satisfy the requirements
                   idx_transect = transect_data(l,:) > SWH;
                    
                    points = lat(lon_out1+l-1,idx_transect);
                    southern_hemi_points = points(points < 0);
                    [~,wid] = size(southern_hemi_points);
                    if wid == 0 % No MIZ
                        distSWH(l,j) = 0;
                    else
                        distSWH(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                    end
                    if printing == 1
                        fprintf(strcat("The width of the SWH MIZ in the ", sector," sector is: %g km\n"), distSWH)
                    end
                end
        
        
        
                % 4.c FSD
                variable = "fsdrad";
                data = data_format(filename,variable,row,lat,lon,dim);
                if plotting == 1
                    figcount = figcount + 1;
                    figure(figcount)
                    [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                        title('FSDrad along transect')
                else
                    [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            
                end
                transect_data = output_data(lon_out1:lon_out2,:); % Take the southern hemisphere
                
        
                for l = 1:lon_out2-lon_out1
                    % Find the cells that satisfy the requirements
                    idx1 = transect_data(l,:) > eps;
                    idx2 = transect_data(l,:) < fsd_max;
                    idx_transect = logical(idx1.*idx2);
                    
                    points = lat(lon_out1+l-1,idx_transect);
                    southern_hemi_points = points(points < 0);
                    [~,wid] = size(southern_hemi_points);
                    if wid == 0 % No MIZ
                        distFSD(l,j) = 0;
                    else
                        distFSD(l,j) = lldistkm([southern_hemi_points(1), l], [southern_hemi_points(end), l]);
                    end
                    if printing == 1
                        fprintf(strcat("The width of the FSD MIZ in the ", sector," sector is: %g km\n"), distFSD)
                    end
                end
                date = update_date(date);
            end
            Data(k,:).SIC = mean(distSIC);
            Data(k,:).SWH = mean(distSWH);
            Data(k,:).FSD = mean(distFSD);
        end
    else % is a transect
        for k = 1:length(months)
        month_init = months(k);
        if month_init < 10
            date = sprintf('%d-0%d-0%d', year, month_init, day);
        else
            date = sprintf('%d-%d-0%d', year, month_init, day);
        end
        if month_init == 1 || month_init == 3 || month_init == 5 || month_init == 7 || month_init == 8 || month_init == 10 ||  month_init == 12
            datapoints = 31;
        elseif month_init == 2
            datapoints = 28;
        else
            datapoints = 30;
        end
        distSIC = zeros(1,datapoints);
        distSWH = zeros(1,datapoints); 
        distFSD = zeros(1,datapoints);
        for j = 1:datapoints
            filename = strcat(filedir,"/history/iceh.",date,".nc");
            % 4. a) SIC
            [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "aice";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
            transect_data = output_data(lon_out,:); % Take the southern hemisphere
            
            
            idx_sic = transect_data > SIC(1);
            idx2 = transect_data < SIC(2);
            idx_transect = logical(idx_sic.*idx2);
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distSIC(j) = 0;
            else
                distSIC(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the SIC MIZ along %g E is %g km\n', sector(2), distSIC)
            end
            
            % 4.b SWH
            [~,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "wave_sig_ht";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
            end
        
            transect_data = output_data(lon_out,:); % Take the southern hemisphere
            

            % Find points where waves are greater than SWH and ice is
            % greater than 15%
            idx_swh = transect_data > SWH;
            idx_transect = logical(idx_sic.*idx_swh); 
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distSWH(j) = 0;
            else
                distSWH(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the SWH MIZ along %g E is %g km\n', sector(2), distSWH)
            end
            % 4.c FSD
            
            [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
            variable = "fsdrad";
            data = data_format(filename,variable,row,lat,lon,dim);
            if plotting == 1
                figcount = figcount + 1;
                figure(figcount)
                [w, a, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
                    title('FSDrad along transect')
            else
                [~, ~, output_data] = map_plot(data,variable,sector,grid,[0,11],plotting);
        
            end
            transect_data = output_data(lon_out,:); % Take the southern hemisphere

            % Find points where FSD < max and SIC > 0.15
            idx1 = transect_data < fsd_max;
            idx2 = transect_data > eps;
            idx_transect = logical(idx1.*idx2.*idx_sic);
            
            points = lat(lon_out,idx_transect);
            southern_hemi_points = points(points < 0);
            [~,wid] = size(southern_hemi_points);
            if wid == 0 % No MIZ
                distFSD(j) = 0;
            else
                distFSD(j) = lldistkm([southern_hemi_points(1), lon_out], [southern_hemi_points(end), lon_out]);
            end
            if printing == 1
                fprintf('The width of the FSD MIZ along %g E is %g km\n', sector(2), distFSD)
            end
            date = update_date(date);
        end
        Data_sic(k).Month = [distSIC;distSWH;distFSD];
    end
    end
end
%% Plotting

if plotting_mean == 1
    if isstring(sector)
        %for i = 1:datapoints
        %    distSIC_ave_daily(m,i) = mean(distSIC(:,i));
        %    distSWH_ave_daily(m,i) = mean(distSWH(:,i));
        %    distFSD_ave_daily(m,i) = mean(distFSD(:,i));
        %end
        % 1a) Boxplot widths
        figcount = figcount + 1;   
        f = figure(figcount);
        t = tiledlayout(1,length(months));
        for i = 1:length(months)
            plot_data = [distSIC_ave_daily;distSWH_ave_daily;distFSD_ave_daily];
            nexttile
            boxplot(plot_data',["SIC","SWH","FSD"])
            ylabel('MIZ width (km)')
            xlabel('MIZ definition')
            ylim([0,ymax])
            text = strcat("MIZ widths of the ", sector, " sector transect over " ,monthName(months(i))," %g");
            title(sprintf(text,year))
            f.Position = [100 100 1500 400];
        end
    
          % 1b) Correlations
        figcount = figcount + 1;
        f = figure(figcount);
        t2 = tiledlayout(1,length(months));
        for i = 1:length(months)
            %plot_data = Data(i).Month;
            nexttile
            x = plot_data(2,:);
            y = plot_data(3,:);
            %scatter(x,y)
            % Get coefficients of a line fit through the data.
            coefficients = polyfit(x, y, 1);
            % Create a new x axis with exactly 1000 points (or whatever you want).
            xFit = linspace(min(x), max(x), 1000);
            % Get the estimated yFit value for each of those 1000 new x locations.
            yFit = polyval(coefficients , xFit);
            % Plot everything.
            plot(x, y, 'b.', 'MarkerSize', 15); % Plot training data.
            hold on; % Set hold on so the next plot does not blow away the one we just drew.
            plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
            grid on;
            ylabel('SWH MIZ width (km)')
            xlabel('FSD MIZ width (km)')
            ylim([0,ymax])
            xlim([0,1000])
            text = strcat(monthName(months(i))," %g in the ", sector, " sector");
            title(sprintf(text,year,sector(2)))
            f.Position = [100 100 1500 400];
        end
        % 1c) Widths over time
        for i = lon_out1:lon_out2
            lengend_text(i-lon_out1+1) = sprintf("%g",i);
        end
        figcount = figcount + 1;
        f = figure(figcount);
        t_width = tiledlayout(1,3);
        nexttile
        plot(1:datapoints,distSIC)
        legend(lengend_text)
        xlabel("Dates")
        ylabel("MIZ width (km)")
        nexttile
        plot(1:datapoints,distSWH)
        xlabel("Dates")
        ylabel("MIZ width (km)")
        nexttile
        plot(1:datapoints,distFSD)
        xlabel("Dates")
        ylabel("MIZ width (km)")
    else
        % 2a) Boxplot widths
        figcount = figcount + 1;
        f = figure(1);
        t = tiledlayout(1,length(months));
        for i = 1:length(months)
            plot_data = Data_sic(i).Month;
            %plot_data = Data(i).Month;
            nexttile
             boxplot(plot_data',["SIC","SWH","FSD"])
%             Data.Month = categorical(Data.Month);
%             load patients
%             binEdges = 25:5:50;
%             bins = {'late 20s','early 30s','late 30s','early 40s','late 40s+'};
%             groupAge = discretize(Age,binEdges,'categorical',bins);
%             boxchart(groupAge,Diastolic)
%             xlabel('Age Group')
%             ylabel('Diastolic Blood Pressure')
% 
%             monthOrder = categorical({'February'});
%             Data.Month = categorical(Data.Month,monthOrder);


%             for i = 1:length(months)
%                 plot_data(i,:) = Data(i).SWH;
%             end
%             boxchart(1:22,plot_data(1,:))
             ylabel('MIZ width (km)')
             xlabel('MIZ definition')
             ylim([0,ymax])
% 
% 
             text = strcat("MIZ widths of the %g E transect over " ,monthName(months(i))," %g");
             title(sprintf(text,sector(2),year))
             f.Position = [100 100 1500 400];
        end
       % 2b) Correlations
        figcount = figcount + 1;
        f = figure(figcount);
        t2 = tiledlayout(1,length(months));
        for i = 1:length(months)
            plot_data = Data(i).Month;
%             nexttile
%             x = plot_data(2,:);
%             y = plot_data(3,:);
%             %scatter(x,y)
%             % Get coefficients of a line fit through the data.
%             coefficients = polyfit(x, y, 1);
%             % Create a new x axis with exactly 1000 points (or whatever you want).
%             xFit = linspace(min(x), max(x), 1000);
%             % Get the estimated yFit value for each of those 1000 new x locations.
%             yFit = polyval(coefficients , xFit);
%             % Plot everything.
%             plot(x, y, 'b.', 'MarkerSize', 15); % Plot training data.
%             hold on; % Set hold on so the next plot does not blow away the one we just drew.
%             plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
%             grid on;
%             ylabel('SWH MIZ width (km)')
%             xlabel('FSD MIZ width (km)')
%             ylim([0,ymax])
%             xlim([0,1000])
%             text = strcat(monthName(months(i))," %g along %g E transect");
%             title(sprintf(text,year,sector(2)))
%             f.Position = [100 100 1500 400];
        end
       % 3a) SIC Boxplots
%         figcount = figcount + 1;
%         f = figure(figcount);
%         t = tiledlayout(1,length(months));
%         for i = 1:length(months)
%             plot_data = Data_sic(i).Month;
%             nexttile
%             boxplot(plot_data',["SIC","SWH and SIC > 0.15","FSD and SIC > 0.15"])
%             
%             %ylabel('Temperature (F)')
%             legend
%             ylabel('MIZ width (km)')
%             xlabel('MIZ definition')
%             ylim([0,ymax])
%             text = strcat("MIZ widths of the %g E transect over " ,monthName(months(i))," %g");
%             title(sprintf(text,sector(2),year))
%             f.Position = [100 100 1500 400];
%         end
%         %3b) SIC Correlations
%          figcount = figcount + 1;
%         f = figure(figcount);
%         t2 = tiledlayout(1,length(months));
%         for i = 1:length(months)
%             plot_data = Data_sic(i).Month;
%             nexttile
%             x = plot_data(2,:);
%             y = plot_data(3,:);
%             %scatter(x,y)
%             % Get coefficients of a line fit through the data.
%             coefficients = polyfit(x, y, 1);
%             % Create a new x axis with exactly 1000 points (or whatever you want).
%             xFit = linspace(min(x), max(x), 1000);
%             % Get the estimated yFit value for each of those 1000 new x locations.
%             yFit = polyval(coefficients , xFit);
%             % Plot everything.
%             plot(x, y, 'b.', 'MarkerSize', 15); % Plot training data.
%             hold on; % Set hold on so the next plot does not blow away the one we just drew.
%             plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
%             grid on;
%             ylabel('SWH and SIC > 0.15 MIZ width (km)')
%             xlabel('FSD and SIC > 0.15 MIZ width (km)')
%             ylim([0,ymax])
%             xlim([0,1000])
%             text = strcat(monthName(months(i))," %g along %g E transect");
%             title(sprintf(text,year,sector(2)))
%             f.Position = [100 100 1500 400];
%        end
    end
end
%% Functions

% function ave_data = aggregate_data(case_name,date,datapoints,variable,dim)
%     if dim == 2
%         for i = 1:datapoints % number of days
%            filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
%            data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
%            data = rearrange_matrix(data,37,dim);
%            data_formatted = [data; data(end,:)];
%            grid = "gx1";
%            mask = ice_mask(case_name,date,grid,SIC);
%            total_data(:,:,i) = data_formatted.*mask;
%            % update date
%            date = update_date(date);
%         end
%         ave_data = mean(total_data,3);
%     else % dim == 3
%         for i = 1:datapoints % number of days
%             filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
%             data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
%             nfsd = ncread(filename, "NFSD");
%             data1D = sum(data,3)/numel(nfsd); % diving by the number of categories
%             data1D = rearrange_matrix(data1D,37,2);
%             total_data(:,:,i) = [data1D; data1D(end,:)];
%             % update date
%             date = update_date(date);
%         end
%         ave_data = mean(total_data,3);
%     end
% end
% 
% function edge = ice_edge(case_name,date,grid)
%     %% Find the ice edge
%     % Find the ice edge and land edge
%     filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
%     latit = 100;
%     variable = "aice";
%     dim = 2;
%     datapoints = 92;
%     [lat,lon,row] = grid_read(grid);
%     aice_data(:,:) = aggregate_data(case_name,date,datapoints,variable,dim);
%     [len, wid] = size(aice_data);
%     for i = 1:len
%         long_ice_edge = aice_data(i,1:latit);
%         pos = find(long_ice_edge < eps);
%         ice_pos(i) = pos(1);
%     end
%     edge = ice_pos;
% end
% 
% 
% function mask = ice_mask(case_name,date,grid,conc)
%     %% Find the ice edge
%     % Find the ice edge and land edge
%     filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
%     latit = 100;
%     variable = "aice";
%     dim = 2;
%     datapoints = 92;
%     [lat,lon,row] = grid_read(grid);
%     filename = strcat("cases/",case_name,"/history/iceh.",date,".nc");
%     aice_data = ncread(filename, variable); % wave_sig_ht, dafsd_wave, fsdrad, peak_period
%     aice_data = rearrange_matrix(aice_data,37,dim);
%     aice_data_formatted = [aice_data; aice_data(end,:)];
%     [len, wid] = size(aice_data_formatted);
%     ice_pos = zeros(len,wid);
%     for i = 1:len
%         long_ice_edge = aice_data_formatted(i,1:latit);
%         pos = find(long_ice_edge > conc);
%         ice_pos(i,pos) = 1;
%     end
%     mask = ice_pos;
% end
% 
% function y = map_creator_miz(filename, plot_title_vec, i, variable, grid, map_type, user, case_name)
% %MAP_CREATOR creates map images of Antarctica given netcdf data files
% %   filename: the directory to the .nc file
% %   plot_title: string containing the title for the plot
% %   i: the date of the data file
% %   variable: the variable in the dataset we want to plot
% %   grid: specify the grid. eg. 'gx1', 'gx3'
% 
% % Load ice shelf data
% addpath packages/bedmap2_toolbox_v4
% 
% shelf = bedmap2_data('icemask');
% shelf = 0.0*(shelf == 1); 
% shelf(shelf==0) = NaN; % 1 iceshelf
% [latshelf,lonshelf] = bedmap2_data('latlon');
% 
% % grid
%     dim = 2;
% if grid == 'gx3'
%     row = 11;
% 
%     ulat = ncread('grid_gx3.nc','ulat');
%     ulon = ncread('grid_gx3.nc','ulon');
% 
%     % converting to degrees
%     lon = rad2deg(ulon);
%     lat = rad2deg(ulat);
% 
%     lat = rearrange_matrix(lat,row,dim);
%     lon = rearrange_matrix(lon,row,dim);
%    
% else
%     row = 37;
%     lat = ncread('grid/global_gx1.bathy.nc','TLAT');
%     lon = ncread('grid/global_gx1.bathy.nc','TLON');
% 
%     lat = rearrange_matrix(lat,row,dim);
%     lon = rearrange_matrix(lon,row,dim);
% 
% 
%     lon = [zeros(1,384);lon];
%     lat = [lat(1,:); lat];
%         
% end
% 
% %filename = strcat('cases/',filename);
% data = ncread(filename, variable);
% [~, ~, n] = size(data);
% 
% for level = 1:n
%     data_1 = data(:,:,level);
% 
%     latitude = [-90,90];
%     longitude = [-180,180];
% 
%     data_1 = rearrange_matrix(data_1,row,dim);
% 
%     % fixing data
%     [m, ~] = size(lon);
%     lon = [lon; lon(end,:) + 360/m];
%     lat = [lat; lat(end,:)];
%     data_1 = [data_1; data_1(end,:)];
%     %% Threshold for the data
%     variable = "fsdrad";
%     fsd_max = 50.0;
%     idx = data_1 > fsd_max; 
%     fsd_miz = data_1;
%     fsd_miz(idx) = 0.0;
%     idx = fsd_miz > eps;
%     fsd_miz(idx) = fsd_max;    
%     %% Mapping
%     color_map = seaicecolormap();
% if map_type == 'cassini'
%     x_origin = -20;
% else
%     x_origin = -90;
% end
%     w = worldmap('world');
%     axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
%     setm(w, 'Origin', [x_origin 0 0]);
%     setm(w, 'maplatlimit', [-90,-30]);
%     setm(w, 'maplonlimit', [-180,180]);
%     setm(w, 'meridianlabel', 'on')
%     setm(w, 'parallellabel', 'off')
%     setm(w, 'mlabellocation', 30);
%     setm(w, 'plabellocation', 10);
%     setm(w, 'mlabelparallel', -45);
%     setm(w, 'grid', 'on');
%     setm(w, 'labelrotation', 'on')
%     pcolorm(lat,lon,fsd_miz)
%     land = shaperead('landareas', 'UseGeoCoords', true);
%     geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
%     pcolorm(latshelf,lonshelf,shelf)  
% 
%  if variable == "fsdrad"
%          plot_variable = "FSD radius ";
%          unit = "metres";
%      elseif  variable == "fsdrad_d"
%          plot_variable = "FSD radius ";
%          unit = "metres";
%      elseif variable == "wave_sig_ht"
%          plot_variable = "Significant wave height ";
%          unit = "metres";
%      elseif variable == "wave_sig_ht_d"
%          plot_variable = "Significant wave height ";
%          unit = "metres";
%      elseif variable == "peak_period"
%          plot_variable = "Peak period ";
%          unit = "s";
%      elseif variable == "peak_period_d"
%          plot_variable = "Peak period ";
%          unit = "s";
%      elseif variable == "aice"
%          plot_variable = "Concentration of ice ";
%          unit = " ";
%      elseif variable == "aice_d"
%          plot_variable = "Concentration of ice ";
%          unit = " ";
%      elseif variable == "mean_wave_dir"
%          plot_variable = "Mean wave direction (rads) ";
%          unit = "radians";
%      elseif variable == "mean_wave_dir_d"
%          plot_variable = "Mean wave direction (rads) ";
%          unit = "radians";
%      else
%          plot_variable = variable;
%  end
%     set(gcf, 'Position',  [100, 100, 1000, 800])
%     set(gcf,'Visible', 'off')
%     fontSize = 20; 
%     plot_title = strcat(plot_variable, plot_title_vec);
%     title(plot_title, 'FontSize', fontSize);
%     a=colorbar;
%     label_c = ylabel(a,unit,'FontSize',16,'Rotation',270);
%     label_c.Position(1) = 4;
%     label_h.Position(2) = 1; % change vertical position of ylabel
%     limit = colorlims(variable);
%     caxis(limit);
%     figname = sprintf('image%d.png', i); 
%     filedir = sprintf('/Users/%s/GitHub/CICE-plotting-tools/frames', user);
%     saveas(gcf,fullfile(filedir, figname));
% end
% end
% 
function [] = plot_map(lat,lon,total_miz,latshelf,lonshelf,shelf,text,i,colourbar)
    %% Mapping
    latitude = [-90,-30];
    longitude = [-180,180];
    addpath functions
    figure(i)
    w = worldmap('world');
        axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
        setm(w, 'Origin', [-90 0 0]);
        setm(w, 'maplatlimit', [-90,-30]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'mlabellocation', 30);
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(latshelf,lonshelf,shelf)
        pcolorm(lat,lon,total_miz)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])  
        colorbar
        if colourbar == 1
            caxis([0,1])
        end
        title(text,'interpreter','latex','FontSize', 18)
        set(gcf, 'Position',  [100, 100, 1000, 800])
end
function name = monthName(num)
    name = month(datetime(1,num,1), 'name');
    name = name{1};
end