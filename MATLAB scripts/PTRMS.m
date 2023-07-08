%%% This is a script for the PTR-MS VOC data %%% 

%% UTC TIME %%
UTC = KORUSAQPTRMSNMHCsDC820160618RA1(:,1); 


%% Latitude/Longitude/Altitude %%

Latitude = KORUSAQHskpingdc820160618RA1(535:20614,3);
Longitude = KORUSAQHskpingdc820160618RA1(535:20614,4);
% Altitude in Meters % 
Altitude = KORUSAQHskpingdc820160618RA1(535:20614,8).*0.3048;

%%% VOCs in pptV %%% 
Acetonitrile_pptV =  KORUSAQPTRMSNMHCsDC820160618RA1(:,2); 
Acetone_Propanal_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,3); 
Toluene_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,4);
Benzene_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,5);
Toluene_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,6);
C8-Aromatics_Benzaldehyde_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,7);
Monoterpenes_pptV = KORUSAQPTRMSNMHCsDC820160618RA1(:,8);

%%% Time %%% 

% UTC is 7 hours ahead of PDT

PDT = UTC - (7 * 60 * 60);

[hrs,min,sec] = seconds2hms(PDT);

start = concatenate(hrs(1),min(1),sec(1));




%%% Temperature and Pressure %%%

% Change from Celsius to Kelvin %
Temperature = KORUSAQHskpingdc820160618RA1(535:20614,19) + 273.15;

% Converts from mb to Pascals % 
Pressure = KORUSAQHskpingdc820160618RA1(535:20614,24).*100;

% Volume is initially in m^3 so the 10^6 is to get it into cm^3% 

Volume = [(10^6)*(1.66 * 10^-12)*(8.314).*Temperature]./Pressure;

%%% ISOPRENE %%% 
Toluene_pptV(Toluene_pptV == -9999) = NaN;
Toluene_pptV(Toluene_pptV < 0) = NaN;
Isoprene_pptV2 = Toluene_pptV;

NumberDensityIsoprene = rdivide(Isoprene_pptV2, Volume); 

% Isoprene Rate Constant with OH %
k_Isoprene = (2.7 * 10^-11)*exp(390./Temperature); 


% Organic Reactivity %
IsopreneOrgReactivity = k_Isoprene.*NumberDensityIsoprene; 

meanIsopreneOR = nanmean(IsopreneOrgReactivity);
sdIsopreneOR = nanstd(IsopreneOrgReactivity); 

%%% Time Series of Isoprene Organic Reactivity %%%

ts1 = timeseries(IsopreneOrgReactivity,1:length(IsopreneOrgReactivity));
ts1.Name = 'Isoprene Organic Reactivity(s^-1)'
ts1.TimeInfo.Units = 'Seconds'
ts1.TimeInfo.StartDate = strcat('18-Jun',start)
ts1.TimeInfo.Format = 'HH:MM:SS';

ts1.Time = ts1.Time - ts1.Time(1)'

plot(ts1,'k'); 
title('6/18/16 Isoprene Organic Reactivity','Interpreter','Tex','FontSize',16);
ylabel('Organic Reactivity s^-1','Interpreter','Tex','FontSize',14);
xlabel('Time (PDT)'); 


%%% TOLUENE %%% 

Toluene_pptV(Toluene_pptV == -9999) = NaN;
Toluene_pptV(Toluene_pptV < 0) = NaN;
Toluene_pptV2 = Toluene_pptV;

%%% 3D Plot %%% 

scatter3(Longitude,Latitude,Altitude,[5],Isopreneppb)
colormap jet
caxis([0,1])
colorbar
ylabel(colorbar, 'ppb','Interpreter', 'tex')

title('Isoprene Concentrations 6/18/2016 ', 'FontSize', 18, 'Interpreter', 'tex')
xlabel('Longitude', 'FontSize', 14)
ylabel('Latitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)


s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 42.5])


%%% Testing Function On Isoprene %%%
IsopreneOrgReac = function_OrgReac(Isoprene_pptV,Temperature,Pressure,(2.7*10^-11)*exp(390./Temperature),PDT,Latitude,Longitude,Altitude);

%% Trying the function on Benzene %% 

BenzeneOrgReac = function_OrgReac(Benzene_pptV, Temperature, Pressure, (2.33 * 10^-12)*exp(-193./Temperature));

%%% Toluene %%%
TolueneOrgReac = function_OrgReac(Toluene_pptV,Temperature,Pressure,1.18*10^-12*exp(338./Temperature),PDT,Latitude,Longitude,Altitude);

%% Trying it with Acetonitrile %%

AcetonitrileOrgReac = function_OrgReac(Acetonitrile_pptV,Temperature,Pressure,(2.3*10^-14),PDT,Latitude,Longitude,Altitude);
