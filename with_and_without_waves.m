clear all
close all
%% 1. Preamble

% Comparing the difference w/ and w/o the WIM module being turned on.
% Also have to do a comparison w/ and w/o the floe size re-initialisation.
% All runs were done over 7 timesteps or 25200 hours (what is this?)

% Case List:
% - casenowaves := waves in ice module turned off
% - casezero := WIM on, with IFD being reset to 0 when ice conc. < puny
% - casenozero := WIM on, commenting out the reset

% Comparing CICE results with and without waves
 total_plot_title = 'Comparison between having WIM on/off';
 cases = [ "8month","1yearnowaves"];%["casezero","casenozero"];%["casenozero", "casenowaves"];
 date = '2005-08-31';
 datapoints = 1;
 grid = 'gx1';
 timestep = 'd'; % '1', 'd', 'm', 'y'
 user = 'noahday'; %a1724548, noahday, Noah
 variable = ["uvel","vvel","hi"]; %["aice","frazil"];%["wave_sig_ht","peak_period","fsdrad"]; 
 pancake = 'off';
 map_type = "eqdcylin"; %cassini, eqdcylin
 % eqaazim = Round Earth
 % eqdcylin = Cylinder, rectangular plot
 rad = 110;
 no_cases = length(cases);
 no_variables = length(variable);
 filedir = [];
 dim = 2;
 
 % sea ice area: 'aicen' or 'aice'
 % horizonal velocity: 'uvel'
 % vertical velocity: 'vvel'
if timestep == '1'
    for i = 1:no_cases
        filedir = [filedir; strcat('cases/', cases(i),'/history/iceh_inst.',date,'.nc')];
    end
else
    for i = 1:no_cases
        filedir = [filedir; strcat('cases/',cases(i),'/history/iceh.',date,'.nc')];
    end
end

% Read the header
ncdisp(filedir(2));

%% 2. Grid and Data

[lat,lon,row] = grid_read(grid);


%% 3. Mapping
    
    t = tiledlayout(no_variables,no_cases+1); % Requires R2019b or later
    idx = 1;
    for j = 1:no_variables
        for i = 1:no_cases+1
           nexttile(idx) 
           if i == 1
             plot_title = variable(j);
           elseif i == no_cases+1
               if j == 1
                plot_title = "Difference";
               end
           else
              plot_title = " ";
           end
           if i ~= no_cases+1
            %if i == 1 
             %   variable(j) = "wave_sig_ht";
            %else
            %    variable = ["wave_sig_ht","peak_period","fsdrad"];
            %end
            %if i == 1
            %    variable = ["wave_sig_ht","fsdrad","aice"];
            %else
            %    variable = ["wave_sig_ht","fsdrad","aice"];
            %end
            data(:,:,i,j) = data_format(filedir(i),variable(j),row,lat,lon,dim);
            if i ==  no_cases
                endline = 1;
            else
                endline = 0;
            end
            
            mapmaker(data(:,:,i,j),lat,lon,plot_title,map_type,variable(j),endline)
           else
            difference_maker(filedir(1), filedir(2),plot_title, variable(j), grid, map_type, user)
           end
           idx = idx+1;
        end
        
    end

        
    
    % Add layout title
    title(t,total_plot_title)
    set(gcf, 'Position',  [100, 100, 1200, 400])
    

    
   