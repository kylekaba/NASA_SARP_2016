
% Inputting altitude, latitude, and longitude for Butene 
% latitude = KORUSAQHskpingdc820160618RA((681:20695),3);
% longitude = KORUSAQHskpingdc820160618RA((681:20695),4);
% altitude = KORUSAQHskpingdc820160618RA((681:20695),5);
% altmeters = floor(altitude/3.2808);

% For Ozone, need to adjust the lengths of the tables
% These are altitudes that correspond to the O3 and NO,NO2,Noy data from
% NCAR
% lat_O3 = KORUSAQHskpingdc820160618RA((512:20011),3);
% long_O3 = KORUSAQHskpingdc820160618RA((512:20011),4);
% alt_O3 =  KORUSAQHskpingdc820160618RA((512:20011),5);
% alt_meters2 = floor(alt_O3/3.2808);

%These altitude values are for the VOCs
lat_VOC = KORUSAQHskpingdc820160618RA1((681:20695),3);
long_VOC = KORUSAQHskpingdc820160618RA1((681:20695),4);
alt_VOC = KORUSAQHskpingdc820160618RA1((681:20695),5).*0.3048;


%Plotting the variables using Claire's code
UTC = korusaqCITHNO3DC820160618RA01S(:,1); 
HNO3 = korusaqCITHNO3DC820160618RA01S(:,2); 


HNO3(HNO3 == -9999) = NaN;
HNO3(HNO3 < 0) = NaN;
HNO3_2 = HNO3; 

scatter3(long_VOC,lat_VOC,alt_VOC,[5],HNO3_2)
colormap jet
caxis([0,1000])
colorbar
ylabel(colorbar, 'HNO_3 (ppt)','Interpreter', 'tex')
% colormap) % 10, 'filled')
% set(gca,'CLim',[0 140]);
% hold on
% plot(day_fraction_o3,ozone_flight)

% labels and format


% This is Ozone's Code 
% title('6/18/2016 Gas Phase HNO_3 Concentrations with Altitude', 'FontSize', '18', 'Interpreter', 'tex')
% xlabel('Latitude', 'FontSize', 14)
% ylabel('Longitude','FontSize', 14)
% zlabel('Altitude (m)','FontSize', 14)

%This is NO_2's code
title('6/18/2016 HNO_3 Concentrations with Altitude', 'FontSize', 18, 'Interpreter', 'tex')
xlabel('Latitude', 'FontSize', 14)
ylabel('Longitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)


s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 42.5])






