% This script will generate a flight map overlayed
% with ozone levels by sensitivity

% %%% Time %%%
% 
% DC8_time = DC8Housekeeping_0618(444:19943,1);
% DC8_flight = table2array(DC8_time);
% day_fraction_dc8 = str2double(DC8_flight)/86400 - 7/24;
% 
% o3_midtime = NOxyO3_0618(:,3);
% o3_flight = table2array(o3_midtime);
% day_fraction_o3 = str2double(o3_flight)/86400 - 7/24;
% % day_double = str2double(cellstr(day_fraction));

%%% Flight Latitude & Longitude %%%

DC8_lat = DC8Housekeeping_0618(444:19943,3);
latitude_dc8 = table2array(DC8_lat);

DC8_long = DC8Housekeeping_0618(444:19943,4);
longitude_dc8 = table2array(DC8_long);

%%% Radar Altitude %%%

DC8_radalt = DC8Housekeeping_0618(444:19943,8);
altitude_dc8_ft = table2array(DC8_radalt);
% convert ft to m
altitude_dc8_m = altitude_dc8_ft*0.3048;

%%% Map Overlay %%%

% import map...???

%%% Ozone Levels %%%

ozone_flight = str2double(NOxyO3_0618{:,7});
% ozone_not = ozone;
ozone_flight(ozone_flight == -999999.9) = NaN;
% ozone_color = linspace(0,140,ozone_flight);
% z = z(:)

%%% Plotting %%%

% create line varying with intensity
% 2D plot
%  t=0:0.1:10*pi;
%  h=colormapline(t.*sin(t),t.*cos(t),[],jet(128));
%  set(h,'linewidth',3,'linestyle','--')

scatter3(longitude_dc8,latitude_dc8,altitude_dc8_m,[.5],ozone_flight)
colormap jet
caxis([0,80])
colorbar
ylabel(colorbar, 'Ozone (ppb)')
% colormap) % 10, 'filled')
% set(gca,'CLim',[0 140]);
% hold on
% plot(day_fraction_o3,ozone_flight)

% labels and format
title('Flight Ozone Concentrations with Altitude')
xlabel('Latitude')
ylabel('Longitude')
zlabel('Altitude (m)')

% datetick('x','HH:MM','keepticks')