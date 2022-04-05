function [dafsd] = read_in_fsd_data(cases,date,datapoints,sector,ssd,timestep)
%READ_IN_DAFSD - Reads in all FSD relevant data and outputs in a struct.
% Ice mask uses a SIC of 1%.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    cases - Matrix of case names
%    date - Initial date (Character)
%    datapoints - Number of days
%    ssd - Whether data is stored ssd or not
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: data_format_sector
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Noah Day
% Work address
% email: noah.day@adelaide.edu.au
% March 2022; Last revision: 29-March-2022

%------------- BEGIN CODE --------------

% INCLUDE DEFAULT STATEMENTS IN CASE INPUTS ARE MISSING


% By default, set the timestepping to daily
    if ~exist('timestep','var')
        timestep = "d";
    end
 % Get FSD info
if ssd == 1
    filename = strcat('/Volumes/NoahDay5TB/cases/',cases(1),'/history/iceh.',date,".nc");
else
    filename = strcat('cases/',cases(1),"/history/iceh.",date,".nc"); 
end
NFSD = ncread(filename,"NFSD");
Nf = numel(NFSD);

lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

initial_date = date;

for j = 1:numel(cases)
    ticker = 1;
    date = initial_date;
    case_name = cases(j);
    for i = 1:datapoints
       % Update the filename
        if ssd == 1
            filename = strcat('/Volumes/NoahDay5TB/cases/',case_name,'/history/iceh.',date,".nc");
        else
            filename = strcat('cases/',case_name,"/history/iceh.",date,".nc"); 
        end
        %aice = data_format_sector(filename,"aice","SH");
        aice = fsd_converter(filename,"afsdn","aice");
        ice_mask = aice > 0.01; % Find all cells with ice greater than 0.01

        % 1) Read in raw data
            temp.latm.raw = data_format_sector(filename,"dafsd_latm",sector); 
            temp.latg.raw = data_format_sector(filename,"dafsd_latg",sector); 
            temp.newi.raw = data_format_sector(filename,"dafsd_newi",sector); 
            temp.weld.raw = data_format_sector(filename,"dafsd_weld",sector); 
            temp.wave.raw = data_format_sector(filename,"dafsd_wave",sector);
            temp.afsd.raw = data_format_sector(filename,"afsd",sector);
    
       % 2) Changes in f(r)
       %    Calculate the average over space    
        for nf = 1:Nf
            % Only taking the cells with SIC > 0.01
            dafsd.latm.raw(:,:,nf,i) = temp.latm.raw(:,:,nf).*ice_mask;
            dafsd.latg.raw(:,:,nf,i) = temp.latg.raw(:,:,nf).*ice_mask;
            dafsd.newi.raw(:,:,nf,i) = temp.newi.raw(:,:,nf).*ice_mask;
            dafsd.weld.raw(:,:,nf,i) = temp.weld.raw(:,:,nf).*ice_mask;
            dafsd.wave.raw(:,:,nf,i) = temp.wave.raw(:,:,nf).*ice_mask; 
            dafsd.afsd.raw(:,:,nf,i) = temp.afsd.raw(:,:,nf).*ice_mask; 
            % Calculate 1/aice mask
            [nxblock, nyblock] = size(aice);
            for nx = 1:nxblock
                for ny = 1:nyblock
                    if isnan(aice(nx,ny))
                        aice_mask(nx,ny) = 0;
                    else
                        aice_mask(nx,ny) = 1/aice(nx,ny);
                    end
                end
            end
            % Take a temporary work variable
            temp.latm.cat(:,:) = dafsd.latm.raw(:,:,nf,i).*aice_mask;
            temp.latg.cat(:,:) = dafsd.latg.raw(:,:,nf,i).*aice_mask;
            temp.newi.cat(:,:) = dafsd.newi.raw(:,:,nf,i).*aice_mask;
            temp.weld.cat(:,:) = dafsd.weld.raw(:,:,nf,i).*aice_mask;
            temp.wave.cat(:,:) = dafsd.wave.raw(:,:,nf,i).*aice_mask;
            temp.afsd.cat(:,:) = dafsd.afsd.raw(:,:,nf,i).*aice_mask;
    
            % Only take the cells where there is more than 1% of SIC
            temp_vec.latm(:) = temp.latm.cat(ice_mask);
            temp_vec.latg(:) = temp.latg.cat(ice_mask);
            temp_vec.newi(:) = temp.newi.cat(ice_mask);
            temp_vec.weld(:) = temp.weld.cat(ice_mask);
            temp_vec.wave(:) = temp.wave.cat(ice_mask);
            temp_vec.afsd(:) = temp.afsd.cat(ice_mask);
                
            % Average over space
            dafsd.latm.ave(nf,i,j) = mean(temp_vec.latm(~isnan(temp_vec.latm)));
            dafsd.latg.ave(nf,i,j) = mean(temp_vec.latg(~isnan(temp_vec.latg)));
            dafsd.newi.ave(nf,i,j) = mean(temp_vec.newi(~isnan(temp_vec.newi)));
            dafsd.weld.ave(nf,i,j) = mean(temp_vec.weld(~isnan(temp_vec.weld)));
            dafsd.wave.ave(nf,i,j) = mean(temp_vec.wave(~isnan(temp_vec.wave)));
            dafsd.afsd.ave(nf,i,j) = mean(temp_vec.afsd(~isnan(temp_vec.afsd)));
 
            clear temp_vec
        end

        % 3) Changes in r_a
            dafsd.latm.ra(i,j) = sum(dafsd.latm.ave(:,i,j).*floe_binwidth');
            dafsd.latg.ra(i,j) = sum(dafsd.latg.ave(:,i,j).*floe_binwidth');
            dafsd.newi.ra(i,j) = sum(dafsd.newi.ave(:,i,j).*floe_binwidth');
            dafsd.weld.ra(i,j) = sum(dafsd.weld.ave(:,i,j).*floe_binwidth');
            dafsd.wave.ra(i,j) = sum(dafsd.wave.ave(:,i,j).*floe_binwidth');
            dafsd.afsd.ra(i,j) = sum(dafsd.afsd.ave(:,i,j).*floe_binwidth');

       % Update date
        date = update_date(date,timestep);
       if mod(i,floor(datapoints/10)) == 0
           clc
           fprintf('%g0%% complete\n',ticker);
           ticker = ticker + 1;
       end
    end    
end

%------------- END OF CODE --------------

%                 % DAFSD converted to changes in r_a
%         temp.ra.latm = fsd_converter(filename,"dafsd_latm","fsdrad",temp.latm.raw,sector);
%         temp.ra.latg = fsd_converter(filename,"dafsd","fsdrad",temp.latg.raw,sector);
%         temp.ra.newi = fsd_converter(filename,"dafsd","fsdrad",temp.newi.raw,sector);
%         temp.ra.weld = fsd_converter(filename,"dafsd","fsdrad",temp.weld.raw,sector);
%         temp.ra.wave = fsd_converter(filename,"dafsd","fsdrad",temp.wave.raw,sector);
%                    % doing it for ra
%             temp_vec.ra.latm(:) = temp.ra.latm(ice_mask);
%             temp_vec.ra.latg(:) = temp.ra.latg(ice_mask);
%             temp_vec.ra.newi(:) = temp.ra.newi(ice_mask);
%             temp_vec.ra.weld(:) = temp.ra.weld(ice_mask);
%             temp_vec.ra.wave(:) = temp.ra.wave(ice_mask);
%             % doing it for ra
%                 dafsd.latm.ra.ave(i) = mean(temp_vec.ra.latm(~isnan(temp_vec.ra.latm)));
%                 dafsd.latg.ra.ave(i) = mean(temp_vec.ra.latg(~isnan(temp_vec.ra.latg)));
%                 dafsd.newi.ra.ave(i) = mean(temp_vec.ra.newi(~isnan(temp_vec.ra.newi)));
%                 dafsd.weld.ra.ave(i) = mean(temp_vec.ra.weld(~isnan(temp_vec.ra.weld)));
%                 dafsd.wave.ra.ave(i) = mean(temp_vec.ra.wave(~isnan(temp_vec.ra.wave)));

        % Change in FSD per day (m/day)
    %     dafsd.latm.ra(:,:,i) = fsd_converter(filename,"dafsd_latm","fsdrad"); 
    %     dafsd.latg.ra(:,:,i) = fsd_converter(filename,"dafsd_latg","fsdrad"); 
    %     dafsd.newi.ra(:,:,i) = fsd_converter(filename,"dafsd_newi","fsdrad"); 
    %     dafsd.weld.ra(:,:,i) = fsd_converter(filename,"dafsd_weld","fsdrad"); 
    %     dafsd.wave.ra(:,:,i) = fsd_converter(filename,"dafsd_wave","fsdrad"); 
    %     
    %     
    %     dafsd_SH.latm.ra(:,:,i) = dafsd.latm.ra(:,:,i)./sector_mask;
    %     dafsd_SH.latg.ra(:,:,i) = dafsd.latg.ra(:,:,i)./sector_mask;
    %     dafsd_SH.newi.ra(:,:,i) = dafsd.newi.ra(:,:,i)./sector_mask;
    %     dafsd_SH.weld.ra(:,:,i) = dafsd.weld.ra(:,:,i)./sector_mask;
    %     dafsd_SH.wave.ra(:,:,i) = dafsd.wave.ra(:,:,i)./sector_mask;
    %     
    %     temp.latm(:,:) = dafsd_SH.latm.ra(:,:,i);
    %     temp.latg(:,:) = dafsd_SH.latg.ra(:,:,i);
    %     temp.newi(:,:) = dafsd_SH.newi.ra(:,:,i);
    %     temp.weld(:,:) = dafsd_SH.weld.ra(:,:,i);
    %     temp.wave(:,:) = dafsd_SH.wave.ra(:,:,i);
    %     
    %     dafsd_SH_ave.latm.ra(i) = mean(temp.latm(~isnan(dafsd_SH.latm.ra(:,:,i))));
    %     dafsd_SH_ave.latg.ra(i) = mean(temp.latg(~isnan(dafsd_SH.latg.ra(:,:,i))));
    %     dafsd_SH_ave.newi.ra(i) = mean(temp.newi(~isnan(dafsd_SH.newi.ra(:,:,i))));
    %     dafsd_SH_ave.weld.ra(i) = mean(temp.weld(~isnan(dafsd_SH.weld.ra(:,:,i))));
    %     dafsd_SH_ave.wave.ra(i) = mean(temp.wave(~isnan(dafsd_SH.wave.ra(:,:,i))));
        
        % Average across the sector
    
        
        