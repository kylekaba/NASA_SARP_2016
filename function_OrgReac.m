function [Organic_Reactivity] = function_OrgReac(pptv, Temperature, Pressure,k)
%time,latitude,longitude,altitude)

% This is a function that will find the Organic Reactivity of a given VOC,
% and also make a time series plot and a 3D plot with respect to the flight
% track. 

%% This is assuming we are in pptv since the 1.66 * 10^-12 represents the amount of moles in 1 trillion air molecules%%
%% Get the volume into cm^3 by multiplying by a factor of 10^6%%

Volume = (1.66*10^-12*8.314.*Temperature*10^6)./Pressure; 

% Making a new vector to store pptv values without -9999 or negative values

A = zeros(length(pptv),1);
A = pptv;

%% Filtering out the sentinel values %% 
for i=1:length(pptv)
if pptv(i,1) == -9999 
    A(i,1) = NaN;
elseif pptv(i,1) < 0
    A(i,1) = NaN; 
    
end
end 
    
NumberDensity = A./Volume; 

Organic_Reactivity = k.*NumberDensity;

% % Now comes in a new function %
% function [hours,minutes,seconds] = seconds2hms(time)
%     hours = floor(time / 3600);
%     time = time - hours*3600;
%     minutes = floor (time/60);
%     seconds = time - minutes*60;
% end
% 
% [hours, minutes, seconds] = seconds2hms(time);
% % And a new function %
% function [concat] = concatenate(hours,minutes,seconds)
%     if seconds < 10;
%         x = strcat('0',num2str(seconds));
%     else 
%         x = num2str(seconds);
%     end
%     
%     if minutes < 10;
%         y = strcat('0',num2str(minutes));
%     else 
%         y = num2str(minutes);
%     end
%     
%     if hours < 10;
%         z = strcat('0',num2str(hours));
%     else 
%         z = num2str(hours);
%     end
%     
%     concat = strcat(z,':',y,':',x);
% end
% 
% start = concatenate(hours(1),minutes(1),seconds(1));
% 
% %Create TimeSeries and Plot
% 
% ts1 = timeseries(Organic_Reactivity,1:length(Organic_Reactivity));
% ts1.Name = 'Organic Reactivity (s^-1)';
% ts1.TimeInfo.Units = 'Seconds';
% ts1.TimeInfo.StartDate = strcat('18-Jun',start);
% ts1.TimeInfo.Format = 'HH:MM:SS';
% 
% ts1.Time = ts1.Time - ts1.Time(1)';
% 
% plot(ts1,'k-'); 
% title('6/18/2016 Organic Reactivity of','Interpreter','Tex','FontSize',16)
% xlabel('Time (PDT)')
% ylabel('Organic Reactivity (s^-1)','Interpreter','Tex')
% 
% 
% %%% Creating a 3D Scatter Plot of the Data %%%
% figure;
% scatter3(longitude,latitude,altitude,[5],Organic_Reactivity)
% colormap jet
% caxis([0,0.5])
% colorbar
% ylabel(colorbar, 'Organic Reactivity','Interpreter', 'Tex')
% 
% 
% 
% title('6/18/2016 Organic Reactivity With Altitude', 'FontSize', 18, 'Interpreter', 'tex')
% xlabel('Longitude', 'FontSize', 14)
% ylabel('Latitude','FontSize', 14)
% zlabel('Altitude (m)','FontSize', 14)
% 
% %%% Superimposing the image onto a map of California %%%
% s=shaperead('CA_counties.shp','UseGeoCoords',true);
% % s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
% hold on
% geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% % geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
% axis([-124 -114 33 42.5])



fprintf('The mean organic reactivity value and standard deviation for this VOC is %g and %g respectively', nanmean(Organic_Reactivity),nanstd(Organic_Reactivity));

end 
















