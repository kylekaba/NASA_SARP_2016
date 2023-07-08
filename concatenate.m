%% Concatenates a time into HH:MM:SS format
function [concat] = concatenate(hours,minutes,seconds)
    if seconds < 10;
        x = strcat('0',num2str(seconds));
    else 
        x = num2str(seconds);
    end
    
    if minutes < 10;
        y = strcat('0',num2str(minutes));
    else 
        y = num2str(minutes);
    end
    
    if hours < 10;
        z = strcat('0',num2str(hours));
    else 
        z = num2str(hours);
    end
    
    concat = strcat(z,':',y,':',x);
end