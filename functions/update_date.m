function new_date = update_date(old_date,timestep)
%UPDATE_DATE - Update the date to the new date.
%
% Syntax:  ew_date = update_date(old_date,timestep)
%
% Inputs:
%    old_date - Character date YYYY-MM-DD
%    timestep - Not essential, "d", "m", "y"
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    date = '2005-07-01';
%    timestep = 'd';
%    new_date = update_date(old_date,timestep)
% ans = '2005-07-02'
%
% Other m-files required: none
% Subfunctions: check_day
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Noah Day
% Email: noah.day@adelaide.edu.au
% March 2022; Last revision: 29-March-2022

%------------- BEGIN CODE --------------
    % By default, set the timestepping to daily
    if ~exist('timestep','var')
        timestep = "d";
    end
    
   
    
    if timestep == "d"
         % Get the initial day, month and year
        current_day = str2num(old_date(end-1:end));
        current_month = str2num(old_date(end-4:end-3));
        current_year = str2num(old_date(end-9:end-6));
        new_day = current_day + 1;
        % Is the new day valid?
        valid_date = check_day(new_day,current_month);
        % Convert to strings
        if current_month < 10
            str_month = sprintf('0%d',current_month);
        else
            str_month = sprintf('%d',current_month);
        end
        if valid_date == 0 % If false update month
            new_day = 1;
            new_month = current_month + 1;
            new_year = current_year;
            if new_month < 10
                str_month = sprintf('0%d',new_month);
            elseif new_month < 13
                str_month = sprintf('%d',new_month);
            else% New year
                new_year = current_year + 1;
                new_month = 1;
                str_month = sprintf('0%d',new_month);
            end
            str_date = strcat('%d-',str_month,'-0%d');
            new_date = sprintf(str_date,new_year,new_day);
        else % If true
            if new_day < 10
                 str_date = strcat('%d-',str_month,'-0%d');
                 new_date = sprintf(str_date,current_year,new_day);
            else
                str_date = strcat('%d-',str_month,'-%d');
                new_date = sprintf(str_date,current_year,new_day);
            end
        end

    elseif timestep == "m"
        % Get the initial month and year
        current_month = str2num(old_date(end-1:end));
        current_year = str2num(old_date(end-6:end-3));
        new_month = current_month + 1;
        % Convert to strings
        if new_month > 12
            % Next year
                 str_date = strcat('%d-','01');
                 new_date = sprintf(str_date,current_year+1);
        else
            % Same year
            if new_month < 10
                str_month = sprintf('0%d',new_month);
            else
                str_month = sprintf('%d',new_month);
            end
             str_date = strcat('%d-',str_month);
             new_date = sprintf(str_date,current_year);
        end
        
    end
end
%------------- END OF CODE --------------

 function dec = check_day(day,month)
    % CHECK_DAY checks if the date is valid
    % dec: is the decision, 1 = valid, 0 = invalid
    
    % NSD: need to add capability for leap years
   
    % What month are we in?
    if  month == 1 || month == 3 || month == 5 || month == 7 || month == 8 || month == 10 || month == 12
        max_day = 31;
    elseif month == 2
        max_day = 28;
    elseif month == 4 || month == 6 || month == 9 || month == 11
        max_day = 30;
    else
        disp('Failed to find month')
        return
    end
    if day < max_day + 1
        dec = 1; 
    else
        dec = 0;
    end
 end
