% Need to write a script that will convert ppm values to molecules/cm^3.
% This will help calculate the organic reactivity.

%Converts the temperature in Celsius to Kelvin. 
Temperature =(KORUSAQHskpingdc820160618RA1(613:20627,19) + 273.15);

%Converts the pressure from mb to Pascals.
Pressure = KORUSAQHskpingdc820160618RA1(613:20627,24)*100;

% CIMS data goes from 61269-81283 UTC.

UTC = korusaqCITPROPENEHNDC820160618RA01S(:,1);


Propene = korusaqCITPROPENEHNDC820160618RA01S(:,2);
Propene(Propene == -9999.0) = NaN; 
Propene(Propene < 0) = NaN;

 

Propene2=Propene;

% To convert to volume we use PV = nRT, V = nRT/P
% P is in Pa, T is in Kelvin, R has units of J·K^?1·mol^?1, n has units of
% moles
% This means V will be in m^3 and multiply this by 1 million too get cm^3.
% 1 trillion air molecules is 1.66 x 10^-12 moles.

%The 10^6 factor is to convert m^3 to cm^3. 
nRT = (1.66 * 10^-12)*(8.314)*(Temperature)*(10^6);

%Volume is now in cm^3
Volume = rdivide(nRT,Pressure);

%Number density is in molecules per cm^3.
NumberDensityPropene = rdivide(Propene2,Volume);

% Defining the rate constant
% Propene's rate constant was found in Pusede (2014)

k_Propene = (4.85*10^-12)*exp(504./Temperature);


% VOCR Organic Reactivity = k_[OH + VOCR][VOCR] 

% Propene's Organic Reactivity
PropeneOrgReactivity = k_Propene.*NumberDensityPropene;

PropeneOrgReacMean = nanmean(PropeneOrgReactivity); 
Propenestd = nanstd(PropeneOrgReactivity);

%%% Latitude Longitude Altitude Data %%%
%%% Going to try and do a 3D plot of Organic Reactivity %%%
lat_O3 = KORUSAQHskpingdc820160618RA((512:20011),3);
long_O3 = KORUSAQHskpingdc820160618RA((512:20011),4);
alt_O3 =  KORUSAQHskpingdc820160618RA((512:20011),5);
alt_meters2 = floor(alt_O3/3.2808);

%Plotting the variables using Claire's code

scatter3(long_O3,lat_O3,alt_meters2,[5],PropeneOrgReactivity(1:19500,1))
colormap jet
caxis([0,0.1])
colorbar
ylabel(colorbar, 'Organic Reactivity (s^-1)','Interpreter','Tex')
% colormap) % 10, 'filled')
% set(gca,'CLim',[0 140]);
% hold on
% plot(day_fraction_o3,ozone_flight)

% labels and format
title('6/18/2016 Propene Organic Reactivity with Altitude', 'FontSize', 18)
xlabel('Latitude', 'FontSize', 14)
ylabel('Longitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)

s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 39])

% datetick('x','HH:MM','keepticks')

%%% Now we will focus on other VOCs %%% 
%%% ETHENE %%%


Ethene(Ethene == -9999) = NaN;
Ethene(Ethene < 0) = NaN;

Ethene2 = Ethene;

%Number Density of Ethene 

NumberDensityEthene = rdivide(Ethene2, Volume);

% Ethene's rate constant

% Ethene's Organic Reactivity
% Need to find literature value for the k-value. 
EtheneOrgReactivity = 

%%% ISOPRENE %%%

% Get rid of all the -9999 and negative values of concentration % 
Isoprene(Isoprene == -9999) = NaN;
Isoprene(Isoprene < 0 ) = NaN;

Isoprene2 = Isoprene;

% Number Density of Isoprene % 
NumberDensityIsoprene = rdivide(Isoprene2, Volume);

%Isoprene rate constants

k_Isoprene = (2.7 * 10^-11)*exp(390./Temperature); 

% VOCR Organic Reactivity = k_[OH + VOCR][VOCR] 

IsopreneOrgReactivity = k_Isoprene.*NumberDensityIsoprene; 

ISOPOrgReacmean = nanmean(IsopreneOrgReactivity)
ISOPstd = nanstd(IsopreneOrgReactivity)

%% Plotting Organic Reactivity With Temperature %%
figure 
subplot(2,2,1);
scatter(Temperature,IsopreneOrgReactivity);
xlabel('Temperature')
ylabel('Organic Reactivity s^{-1','Interpreter','tex','FontSize',16)
title('Isoprene','FontSize',18)

subplot(2,2,2);
scatter(Temperature,PropeneOrgReactivity);
xlabel('Temperature')
ylabel('Organic Reactivity s^{-1}','Interpreter','tex','FontSize',16)
title('Propene','FontSize',18) 










