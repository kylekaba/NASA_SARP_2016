% Getting rid of all the values that are bad
NO_NEW = NO;
NO(NO ==  -999999.9) = NaN;

NO2_2 = NO2;
NO2(NO2 == -999999.9) = NaN;

NOy_2 = NOy;
NOy(NOy == -999999.9) = NaN;

O3_2 = O3;
O3(O3 == -999999.9) = NaN;

% Converting UTC to PST
UTC = MidTime;

PDT = UTC - (7 * 3600);

[hrs, min, sec] = seconds2hms(PDT);

start = concatenate(hrs(1),min(1),sec(1));

% Creating the Time Series Plots
 
ts1 = timeseries(NO_NEW, 1:length(NO_NEW)); 
ts1.Name = ('NO Concentration')
ts1.TimeInfo.Units = ('Seconds')
ts1.TimeInfo.StartDate = strcat('17-Jun',start)
ts1.TimeInfo.Format = 'HH:MM:SS';

ts1.Time = ts1.Time - ts1.Time(1)'

ts2 = timeseries(NO2_2, 1:length(NO2_2));
ts2.Name = ('NO_2 Concentration')
ts2.TimeInfo.Units = ('Seconds')
ts2.TimeInfo.StartDate = strcat('17-Jun', start)
ts2.TimeInfo.Format = 'HH:MM:SS';

ts2.Time = ts2.Time - ts2.Time(1)'

ts3 = timeseries(NOy_2, 1:length(NOy_2));
ts3.Name = ('NO_y Concentration')
ts3.TimeInfo.Units = ('Seconds')
ts3.TimeInfo.StartDate = strcat('17-Jun',start)
ts3.TimeInfo.Format = 'HH:MM:SS';

ts3.Time = ts3.Time - ts3.Time(1)'

ts4 = timeseries(O3_2, 1:length(O3_2));
ts4.Name = ('O_3 Concentration') 
ts4.TimeInfo.Units = ('Seconds')
ts4.TimeInfo.StartDate = strcat('17-Jun',start)
ts4.TimeInfo.Format = 'HH:MM:SS';

ts4.Time = ts4.Time - ts4.Time(1)'

figure
plot(ts1,'r')
hold on

plot(ts2,'b')

hold on

plot(ts3, 'g')

hold on

yyaxis right 
ylabel 'ppb (Volume)'

plot(ts4, 'k')


xlabel('Time (PDT)')
yyaxis left
ylabel 'ppt (Volume)'
yyaxis right 
ylabel 'ppb (Volume)'
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')


legend('NO_2', 'NO', 'NO_y', 'O_3')
title('6/17/2016 NO, NO_2, and NO_y (pptv) O_3 concentration (ppbv)')



