
% Inputting altitude, latitude, and longitude for Butene 
% latitude = KORUSAQHskpingdc820160618RA((681:20695),3);
% longitude = KORUSAQHskpingdc820160618RA((681:20695),4);
% altitude = KORUSAQHskpingdc820160618RA((681:20695),5);
% altmeters = floor(altitude/3.2808);

% For Ozone, need to adjust the lengths of the tables 
lat_O3 = KORUSAQHskpingdc820160618RA((512:20011),3);
long_O3 = KORUSAQHskpingdc820160618RA((512:20011),4);
alt_O3 =  KORUSAQHskpingdc820160618RA((512:20011),5);
alt_meters2 = floor(alt_O3/3.2808);


%Plotting the variables using Claire's code

scatter3(long_O3,lat_O3,alt_meters2,[5],Ozone)
colormap jet
caxis([0,80])
colorbar
ylabel(colorbar, 'Ozone (ppb)')
% colormap) % 10, 'filled')
% set(gca,'CLim',[0 140]);
% hold on
% plot(day_fraction_o3,ozone_flight)

% labels and format
title('6/18/2016 Ozone Concentrations with Altitude', 'FontSize', 18)
xlabel('Latitude', 'FontSize', 14)
ylabel('Longitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)

% datetick('x','HH:MM','keepticks')

