% Plotting data from KORUSAQ
% Need to rescale O3 is in ppb not ppt
% plot(MidTime,NO2pptv,'r',MidTime,O3ppbv,'b',MidTime,NOpptv,'g',MidTime,Noy_ppt,'k')


Midtime = KORUSAQNOxyO3DC820160618RB1(:,1);
NO2pptv = KORUSAQNOxyO3DC820160618RB1(:,4);
Noy_ppt = KORUSAQNOxyO3DC820160618RB1(:,3);
NOpptv = KORUSAQNOxyO3DC820160618RB1(:,2);
O3ppbv = KORUSAQNOxyO3DC820160618RB1(:,5);

% Attempting to Filter out the sentinels 

NO2pptv_2=NO2pptv;
NO2pptv(NO2pptv == -999999.9) = NaN;

NOpptv_2=NOpptv;
NOpptv_2(NOpptv == -999999.9) = NaN;

Noy_ppt_2=Noy_ppt;
Noy_ppt_2(Noy_ppt == -999999.9) = NaN;

O3ppbv_2=O3ppbv;
O3ppbv(O3ppbv == -999999.9) = NaN;


%%% Temperature and Pressure %%% 
%%% Converting to Kelvin and Pa %%% 
Temperature = KORUSAQHskpingdc820160618RA1(444:19943,19) + 273.15;
Pressure = KORUSAQHskpingdc820160618RA1(444:19943,24).*100; 


%% Latitude, Longitude, and Altitude %% 
Latitude = KORUSAQHskpingdc820160618RA1(444:19943,3);
Longitude = KORUSAQHskpingdc820160618RA1(444:19943,4);
Altitude = KORUSAQHskpingdc820160618RA1(444:19943,8).*0.3048;

%%% Calculating the air density for a given function %%%

AirDensity = fun_air(Temperature,Pressure); 

%%% Calculating NO2 and NO number density in molecules/cm^3 %%%
% Mixing ratio of NO2 is on the order of 2*10^-6 percent so multiply by 2*10^-8%

NO2_NumberDensity = (2*10^-8).*(AirDensity);


%%% Calculating the temperature dependent rate constant of OH with NO2 %%%

NO2_OHRateConstant = kohno2new(Temperature,AirDensity);

NO2_OHRateConstantMean = mean(NO2_OHRateConstant);

fprintf('The mean rate value is %g,', mean(NO2_OHRateConstant)) 

% Rate Constant %


% Testing to see if I get the same value both ways, with my function
% opposed to Sally's 2 functions

% Need to multiply the above Number Density with each rate constant %

NO2_OrganicReactivity_SecondTry = NO2_NumberDensity.*NO2_OHRateConstant;
nanmean(NO2_OrganicReactivity_SecondTry) 

%% Using self made function to find NO2's organic reactivity with OH %% 
NO2_OrganicReactivity = function_OrgReac(NO2pptv,Temperature,Pressure,NO2_OHRateConstantMean,PDT,Latitude,Longitude,Altitude); 

%% Now with NO, using Sally's function %%%

NO_RateConstant = knooh(Temperature,AirDensity); 

NO_OrganicReactivity = function_OrgReac(NOpptv_2, Temperature, Pressure, NO_RateConstant, PDT, Latitude,Longitude,Altitude);


%% Combining the two %%
NOxtotal = NO2_OrganicReactivity + NO_OrganicReactivity; 






% UTC is 7 hours ahead of PDT


PDT = Midtime - (7 * 60 * 60);

[hrs,min,sec] = seconds2hms(PDT);

start = concatenate(hrs(1),min(1),sec(1));


%Create time series and plot 

% plot(MidTime,NO2pptv_2,'r',MidTime,NOpptv_2,'g',MidTime,Noy_ppt_2,'b')
% Keeping this here in the event I need to get rid of the previous thing I
% did 

%Create TimeSeries and Plot

ts1 = timeseries(NO2pptv_2,1:length(NO2pptv_2));
ts1.Name = 'NO_2 Concentration (pptv)'
ts1.TimeInfo.Units = 'Seconds'
ts1.TimeInfo.StartDate = strcat('18-Jun',start)
ts1.TimeInfo.Format = 'HH:MM:SS';

ts1.Time = ts1.Time - ts1.Time(1)'

ts2 = timeseries(NOpptv_2,1:length(NOpptv_2));
ts2.Name = 'NO Concentration (pptv)'
ts2.TimeINfo.Units = 'Seconds'
ts2.TimeInfo.StartDate = strcat('18-Jun', start)
ts2.TimeInfo.Format = 'HH:MM:SS';

ts2.Time = ts2.Time - ts2.Time(1)'

ts3 = timeseries(Noy_ppt_2,1:length(Noy_ppt_2));
ts3.Name = 'No_y Concentration (pptv)'
ts3.TimeInfo.Units = 'Seconds'
ts3.TimeInfo.StartDate = strcat('18-Jun', start)
ts3.TimeInfo.Format = 'HH:MM:SS';

ts3.Time = ts3.Time - ts3.Time(1)'


% Plotting ozone separately from everything else 

ts4 = timeseries(O3ppbv,1:length(O3ppbv))
ts4.Name = 'O3 Concentration (ppbv)'
ts4.TimeInfo.Units = 'Seconds'
ts4.TimeInfo.StartDate = strcat('18-Jun', start)
ts4.TimeInfo.Format = 'HH:MM:SS';

ts4.Time = ts4.Time - ts4.Time(1)'

% figure
% plot(ts4, 'k')




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

legend('NO_2', 'NO', 'NO_y', 'O_3')
title('6/18/2016 NO, NO_2, and NO_y (pptv) O_3 concentration (ppbv)')


% Keeping this here in the event I need it.
% xlabel('MidPoint Time in UT seconds')
% ylabel('Parts Per Trillion')
% title('Nitrogen Oxide, Nitrogen Dioxide, and Total Reactive Nitrogen Mixing Ratio')


%% Looking at NOx Organic Reactivity %%

scatter3(Longitude,Latitude,Altitude,[5],NOxtotal)
colormap jet
caxis([0,1])
colorbar
ylabel(colorbar, 'Organic Reactivity s^-1','Interpreter', 'tex')
xlabel('Longitude', 'FontSize', 14)
ylabel('Latitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)
title('Total NOx Reactivity With OH 6/18/2016','FontSize',24) 
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 39.5])

% Time Series For NOx Organic Reactivity % 

OHReac = korusaqmrg01dc8merge20160618RH1(745:20244,1);
OHReac(OHReac == -999999) = NaN;

%% Subtracting NOx reactivity from OH reactivity %% 
OHandNOxDifference = OHReac - NOxtotal; 

%% NOx in ppb %%

NOxppb = (NO2pptv_2 + NOpptv_2).*0.001; 
