function new_date = update_date(old_date)
    % Update the date
    % old_date: character array YYYY-MM-DD
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
        if new_month < 10
            str_month = sprintf('0%d',new_month);
        end
        str_date = strcat('%d-',str_month,'-0%d');
        new_date = sprintf(str_date,current_year,new_day);
    else % If true
        if new_day < 10
             str_date = strcat('%d-',str_month,'-0%d');
             new_date = sprintf(str_date,current_year,new_day);
        else
            str_date = strcat('%d-',str_month,'-%d');
            new_date = sprintf(str_date,current_year,new_day);
        end
    end
 end
 
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
        return
    end
    if day < max_day + 1
        dec = 1; 
    else
        dec = 0;
    end
 end