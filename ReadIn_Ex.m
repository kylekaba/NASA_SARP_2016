%% Read-in Data Example
% Data was read in using the "Import Data" button under the Home tab

UTC = VarName1;
CO2raw = VarName2;

% UTC is 7 hours ahead of PDT

PDT = UTC - (7 * 60 * 60);

[hrs,min,sec] = seconds2hms(PDT);
%start = strcat(num2str(hrs(1)),':',num2str(min(1)),':',num2str(sec(1))); 
start = concatenate(hrs(1),min(1),sec(1));

%% Filter out sentinels

CO2 = CO2raw;
CO2(CO2raw == -9999) = NaN;

CO3 = CO2 + 50;

%% Create Timeseries and plot

ts1 = timeseries(CO2,1:length(CO2));

ts1.Name = 'CO_2 Concentration (ppmv)';
ts1.TimeInfo.Units = 'Seconds';
ts1.TimeInfo.StartDate = strcat('17-Jun',start);
ts1.TimeInfo.Format = 'HH:MM:SS';

ts1.Time = ts1.Time - ts1.Time(1);

ts2 = timeseries(CO3,1:length(CO3));

ts2.Name = 'CO_2 Concentration (ppmv)';
ts2.TimeInfo.Units = 'Seconds';
ts2.TimeInfo.StartDate = strcat('17-Jun',start);
ts2.TimeInfo.Format = 'HH:MM:SS';

ts2.Time = ts2.Time - ts2.Time(1);

figure
plot(ts1)
hold on
plot(ts2)
xlabel('Time (PDT)')
ylabel('CO_2 Concentration (ppmv)')
title('6/17/2016 CO_2 Concentration')






