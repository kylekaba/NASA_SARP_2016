%% Converts seconds after midnight into hours, minutes, seconds

function [hours,minutes,seconds] = seconds2hms(t)
    hours = floor(t / 3600);
    t = t - hours*3600;
    minutes = floor (t/60);
    seconds = t - minutes*60;
end

