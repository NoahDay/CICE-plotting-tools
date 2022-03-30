function [Ticks,TickLabels] = log_ticks(maximum,num_of_ticks,setting)
%LOG_TICKS - Creates colorbar tikz and ticklabels for log transformed data.
%
% Syntax:  [Ticks,TickLabels] = log_ticks(max,num_of_ticks,setting)
%
% Inputs:
%    max - Maximum value of the data
%    num_of_ticks - Number of tick markers wanted
%    setting - Setting for the range of the colorbar
%
% Outputs:
%    Ticks - Tick markers
%    Ticklabels - Associated tick labels
%
% Example: 
%    Data = rand(321,384);
%    maximum = max(Data,[],'all')
%    num_of_ticks = 10;
%    setting = 'balance' % Have 5 data points on either side of 0
%    [Ticks,TickLabels] = log_ticks(max,num_of_ticks,setting)
%    [figs.wave,a] = map_plot(Data,"dafsd_wave",sector);  
%    a.Ticks = Ticks;
%    a.TickLabels = TickLabels;
%    caxis([-log(maximum) log(maximum)]);
%
% Other m-files required: Can be used with map_plot.
% Subfunctions: none
% MAT-files required: none
%
% See also: fstd_drivers.m for a use case.

% Author: Noah Day
% email: noah.day@adelaide.edu.au
% March 2022; Last revision: 25-March-2022

%------------- BEGIN CODE --------------
    % get the minimum and maximum value of A
    if setting == "balance"
        % Let the middle of the cbar be 0
        c1 = 0;
        c2 = maximum;
        % preallocate Ticks and TickLabels
        Ticks      = zeros(1,num_of_ticks);
        TickLabels = zeros(1,num_of_ticks);
        % distribute Ticks and TickLabels
        for n = 1:num_of_ticks/2
            Ticks(num_of_ticks/2+1-n)      = -log(c2)*n/(num_of_ticks/2);
            TickLabels(num_of_ticks/2+1-n) = -round(c2)*n/(num_of_ticks/2);
        end
        for n = 1:num_of_ticks/2
            Ticks(n+num_of_ticks/2)      = log(c2)*n/(num_of_ticks/2);
            TickLabels(n+num_of_ticks/2) = round(c2)*n/(num_of_ticks/2);
        end
        
        Ticks = [Ticks(1:num_of_ticks/2), 0, Ticks(num_of_ticks/2+1:num_of_ticks)];
        TickLabels = [TickLabels(1:num_of_ticks/2), 0, TickLabels(num_of_ticks/2+1:num_of_ticks)];
    elseif setting == "10"
        % Log 10
         % Let the middle of the cbar be 0
        c1 = 0;
        c2 = maximum;
        % preallocate Ticks and TickLabels
        Ticks      = zeros(1,num_of_ticks);
        TickLabels = zeros(1,num_of_ticks);
        % distribute Ticks and TickLabels
        for n = 1:num_of_ticks
            Ticks(n)      = log10(c2)*n/(num_of_ticks);
            TickLabels(n) = round(c2)*n/(num_of_ticks);
        end
        %Ticks = [Ticks(1:num_of_ticks/2), 0, Ticks(num_of_ticks/2+1:num_of_ticks)];
        %TickLabels = [TickLabels(1:num_of_ticks/2), 0, TickLabels(num_of_ticks/2+1:num_of_ticks)];
    end

end
%------------- END CODE --------------
