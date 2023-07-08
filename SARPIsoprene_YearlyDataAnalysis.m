%%% 2016 SARP WAS DATA %%%
Latitude2016 = SARP2016WASDataFlight2(:,11);
Longitude2016 = SARP2016WASDataFlight2(:,12);
Altitude2016 = SARP2016WASDataFlight2(:,14).*0.3048; 

numofdecimalpoints = 1; 

roundedLatitude2016 = round(Latitude2016,numofdecimalpoints); 
roundedLongitude2016 = round(Longitude2016,numofdecimalpoints); 

%%% Now looking at Isoprene levels from that flight %%% 
Isoprene2016 = SARP2016WASDataFlight2(:,92); 
Isoprene2016(Isoprene2016 == -888) = NaN;


%%% Now looking at Alpha Pinene levels from that flight %%%
AlphaPinene2016 = SARP2016WASDataFlight2(:,122); 


%% Isoprene Spatial Variability in 2016 %%
%%% For Loop for Isoprene Plot 

for i = 1:length(Latitude2016); % all points in 2009 scatterp lot
    if Latitude2016(i) >= 35.1 && Latitude2016(i) <= 37.3 && Longitude2016(i) >= -119.9 && Longitude2016(i) <= -119.2;

        subsetLatitude2016(tick) = roundedLatitude2016(i);
        subsetLongitude2016(tick) = roundedLongitude2016(i);
        subsetIsoprene2016(tick) = Isoprene2016(i);
        subsetAltitude2016(tick) = Altitude2016(i); 
        subsetAlphaPinene2016(tick) = AlphaPinene2016(i); 
        tick = tick + 1;
    end
end 
        


scatter3(subsetLongitude2016,subsetLatitude2016,subsetAltitude2016,[80],subsetIsoprene2016,'filled')
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',20)
xlabel('Longitude', 'FontSize', 20)
ylabel('Latitude','FontSize', 20)
zlabel('Altitude (m)','FontSize', 20)
title('Isoprene Levels 6/18/2016','FontSize',30)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])


%% Looking at the SARP 2015 Data File %%

Latitude2015 = SARP2015June24Flight(:,7);
Longitude2015 = SARP2015June24Flight(:,8);
Altitude2015 = SARP2015June24Flight(:,12);

roundedLatitude2015 = round(Latitude2015,numofdecimalpoints); 
roundedLongitude2015 = round(Longitude2015,numofdecimalpoints); 

%% Isoprene in pptv %% 
Isoprene2015 = SARP2015June24Flight(:,74);
Isoprene2015(Isoprene2015 == -888) = NaN;
Isoprene2015(Isoprene2015 < 0) = NaN; 
Isoprene2015_2 = Isoprene2015; 

%%% Looking at Alpha-Pinene %%% 
AlphaPinene2015 =  SARP2015June24Flight(:,96);
AlphaPinene2015(AlphaPinene2015 == -888) = NaN; 

%%% For Loop for 2015 %%% tick = 1;
for i = 1:length(Latitude2015); % all points in 2009 scatterp lot
    if Latitude2015(i) >= 35.1 && Latitude2015(i) <= 37.3 && Longitude2015(i) >= -119.9 && Longitude2015(i) <= -119.2;

        subsetLatitude2015(tick) = roundedLatitude2015(i);
        subsetLongitude2015(tick) = roundedLongitude2015(i);
        subsetIsoprene2015(tick) = Isoprene2015_2(i);
        subsetAltitude2015(tick) = Altitude2015(i); 
        subsetAlphaPinene2015(tick) = AlphaPinene2015(i); 
        tick = tick + 1;
    end
end 
        

%%%% Now Looking at SARP 2014 Data %%%%

%% SARP Flight on 6/24/2014 %% 
% Isoprene6_24_14 = June242014Flight1(:,28);
% Isoprene6_24_14(Isoprene6_24_14 == -888) = NaN;
% Isoprene6_24_14(Isoprene6_24_14 < 0) = NaN;
% Isoprene6_24_2014_2 = Isoprene6_24_14; 
% 
% Latitude6_24_2014 = June242014Flight1(:,9);
% Longitude6_24_2014 = June242014Flight1(:,10);
% Altitude6_24_2014 = June242014Flight1(:,11).*0.3048;
% 
% %% Rounded stuff %% 
% roundedLatitude6_24_2014 = round(Latitude6_24_2014,numofdecimalpoints); 
% roundedLongitude6_24_2014 = round(Longitude6_24_2014,numofdecimalpoints);


%%% 6/25/2014 %%%
Isoprene6_25_14 = June252014(:,28);
Isoprene6_25_14(Isoprene6_25_14 == -888) = NaN;
Isoprene6_25_14(Isoprene6_25_14 < 0) = NaN;
Isoprene6_25_14_2 = Isoprene6_25_14; 

Latitude6_25_2014 = June252014(:,9);
Longitude6_25_2014 = June252014(:,10);
Altitude6_25_2014 = June252014(:,11).*0.3048; 



%%% Alpha Pinene Levels %%%
AlphaPinene6_25_14 = June252014(:,50);
AlphaPinene6_25_14(AlphaPinene6_25_14 == -888) = NaN; 

%% rounded Lat/Longs %%
roundedLatitude6_25_2014 = round(Latitude6_25_2014,numofdecimalpoints);
roundedLongitude6_25_2014 = round(Longitude6_25_2014,numofdecimalpoints);

%%% For loop for 2014 Flights %%%
tick = 1;
for i = 1:length(Latitude6_25_2014); % all points in 2009 scatterp lot
    if Latitude6_25_2014(i) >= 35.1 && Latitude6_25_2014(i) <= 37.3 && Longitude6_25_2014(i) >= -119.9 && Longitude6_25_2014(i) <= -119.2;

        subsetLatitude6_25_2014(tick) = roundedLatitude6_25_2014(i);
        subsetLongitude6_25_2014(tick) = roundedLongitude6_25_2014(i);
        subsetIsoprene6_25_2014(tick) = Isoprene6_25_14_2(i);
        subsetAltitude6_25_2014(tick) = Altitude6_25_2014(i); 
        subsetAlphaPinene6_25_2014(tick) = AlphaPinene6_25_14(i); 
        tick = tick + 1;
    end
end 
        

%%%% Now Looking at SARP 2013 Data %%%%
% 6/18/2013 Flight % 

% Latitude6_18_2013 = SARP6182013Flight(:,9);
% Longitude6_18_2013 = SARP6182013Flight(:,10);
% Altitude6_18_2013 = SARP6182013Flight(:,11).*0.3048;
% 
% roundedLatitude6_18_2013 = round(Latitude6_18_2013,numofdecimalpoints);
% roundedLongitude6_18_2013 = round(Longitude6_18_2013,numofdecimalpoints); 

% 6/19/2013 Flight %
Latitude6_19_2013 = SARP6192013Flight(:,9);
Longitude6_19_2013 = SARP6192013Flight(:,10);
Altitude6_19_2013 = SARP6192013Flight(:,11).*0.3048;

roundedLatitude6_19_2013 = round(Latitude6_19_2013,numofdecimalpoints);
roundedLongitude6_19_2013 = round(Longitude6_19_2013,numofdecimalpoints);

%% Now looking at Isoprene levels from the flight %% 
Isoprene6_19_2013 = SARP6192013Flight(:,63);
Isoprene6_19_2013(Isoprene6_19_2013 == -888) = NaN; 


tick = 1;
for i = 1:length(Latitude6_19_2013); % all points in 2009 scatterp lot
    if Latitude6_19_2013(i) >= 35.1 && Latitude6_19_2013(i) <= 37.3 && Longitude6_19_2013(i) >= -119.9 && Longitude6_19_2013(i) <= -119.2;

        subsetLatitude6_19_2013(tick) = Latitude6_19_2013(i);
        subsetLongitude6_19_2013(tick) = Longitude6_19_2013(i);
        subsetIsoprene6_19_2013(tick) = Isoprene6_19_2013(i);
        subsetAltitude6_19_2013(tick) = Altitude6_19_2013(i); 
        tick = tick + 1;
    end
end 
        

% Scatter Plot of the Isoprene % 
% scatter3(roundedLongitude6_19_2013,roundedLatitude6_19_2013,Altitude6_19_2013,[40],Isoprene6_19_2013,'filled')
% s=shaperead('CA_counties.shp','UseGeoCoords',true);
% % s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
% hold on
% geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% % geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
% axis([-119.9 -119.2 35.1 37.3])




%%% Now Looking at SARP 2012 Flight Data %%%
% 6/25/2012 Flight %
% 
% Latitude6_25_2012 = SARP2012625F1(:,13);
% Longitude6_25_2012 = SARP2012625F1(:,14);
% Altitude6_25_2012 = SARP2012625F1(:,15);
% 
% roundedLatitude6_25_2012 = round(Latitude6_25_2012,numofdecimalpoints);
% roundedLongitude6_25_2012 = round(Longitude6_25_2012,numofdecimalpoints);
% 
% 
% 
% scatter3(Longitude6_25_2012,Latitude6_25_2012,Altitude6_25_2012)
% s=shaperead('CA_counties.shp','UseGeoCoords',true);
% % s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
% hold on
% geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% % geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
% axis([-124 -114 33 42.5])


%%% Now looking at SARP 2011 Flight Data %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Looking at SARP 2010 Data %%% 
Latitude7_1_2010 = SARP201071(:,10);
Longitude7_1_2010 = SARP201071(:,11);
Altitude7_1_2010 = SARP201071(:,12).*0.3048;

roundedLatitude7_1_2010 = round(Latitude7_1_2010,numofdecimalpoints);
roundedLongitude7_1_2010 = round(Longitude7_1_2010,numofdecimalpoints);

%% Isoprene Levels From 2010 Flight %% 

Isoprene7_1_2010 = SARP201071(:,66);
Isoprene7_1_2010(Isoprene7_10_2010 == -888) = NaN; 
Isoprene7_1_2010_2 = Isoprene7_1_2010;


%%% Alpha-Pinene Levels From 2010 %%% 
AlphaPinene7_1_2010 = SARP201071(:,26);

%%% For loop that will filter out only the values defined in this grid. 

tick = 1;
for i = 1:length(Latitude7_1_2010); % all points in 2009 scatterp lot
    if Latitude7_1_2010(i) >= 35.1 && Latitude7_1_2010(i) <= 37.3 && Longitude7_1_2010(i) >= -119.9 && Longitude7_1_2010(i) <= -119.2;

        subsetLatitude7_1_2010(tick) = roundedLatitude7_1_2010(i);
        subsetLongitude7_1_2010(tick) = roundedLongitude7_1_2010(i);
        subsetIsoprene7_1_2010(tick) = Isoprene7_1_2010_2(i);
        subsetAltitude7_1_2010(tick) = Altitude7_1_2010(i); 
        subsetAlphaPinene7_1_2010(tick) = AlphaPinene7_1_2010(i); 
        tick = tick + 1;
    end
end 
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Looking at SARP 2009 Data %%% 
Latitude7_23_2009 = SARP7232009F2(:,13);
Longitude7_23_2009 = SARP7232009F2(:,16);
Altitude7_23_2009 = SARP7232009F2(:,17).*0.3048;

%%% Rounded 2009 Lat/Long %%%
roundedLatitude7_23_2009 = round(Latitude7_23_2009,numofdecimalpoints);
roundedLongitude7_23_2009 = round(Longitude7_23_2009,numofdecimalpoints); 



%% 2009 Isoprene Levels %%
Isoprene7_23_2009 =  SARP7232009F2(:,78);
Isoprene7_23_2009(Isoprene7_23_2009 == -888) = NaN; 


%%% Alpha Pinene levels in 2009 %% 
AlphaPinene7_23_2009 = SARP7232009F2(:,39);


%%% For Loop Let's See %%%

tick = 1;
% subsetLatitude7_23_2009 = zeros(200,1);
%  subsetLongitude7_23_2009 = zeros(200,1);
%  subsetIsoprene7_23_2009 = zeros(200,1); 

for i = 1:length(Latitude7_23_2009); % all points in 2009 scatterp lot
    if Latitude7_23_2009(i) >= 35.1 && Latitude7_23_2009(i) <= 37.3 && Longitude7_23_2009(i) >= -119.9 && Longitude7_23_2009(i) <= -119.2;

        subsetLatitude7_23_2009(tick) = roundedLatitude7_23_2009(i);
        subsetLongitude7_23_2009(tick) = roundedLongitude7_23_2009(i);
        subsetIsoprene7_23_2009(tick) = Isoprene7_23_2009(i);
        subsetAltitude7_23_2009(tick) = Altitude7_23_2009(i); 
        subsetAlphaPinene7_23_2009(tick) = AlphaPinene7_23_2009(i); 
        tick = tick + 1;
    end
end 
        

%% MINTERSECT TIME %%%

mintersect(roundedLatitude2016,roundedLatitude2015,roundedLatitude6_25_2014,roundedLatitude6_19_2013,roundedLatitude7_1_2010,roundedLatitude7_23_2009)

mintersect(roundedLongitude2016,roundedLongitude2015,roundedLongitude6_25_2014,roundedLongitude6_19_2013,roundedLongitude7_1_2010,roundedLongitude7_23_2009) 

%%%%%%%%%%%%%%%%%%%%%% All the Spatial Variability %%%%%

%%%% Changing from ppt to ppb %%%

subsetIsoprene7_23_2009ppb = subsetIsoprene7_23_2009.*0.001; 
subsetIsoprene7_1_2010ppb = subsetIsoprene7_1_2010.*0.001;
subsetIsoprene6_25_2014ppb = subsetIsoprene6_25_2014.*0.001;
subsetIsoprene2015ppb = subsetIsoprene2015.*0.001;
subsetIsoprene2016ppb = subsetIsoprene2016.*0.001;


%%%% Averages, Standard Deviations, %%%%

meanIsoprene2009 = nanmean(subsetIsoprene7_23_2009ppb);
meanIsoprene2010 = nanmean(subsetIsoprene7_1_2010ppb);
meanIsoprene2014 = nanmean(subsetIsoprene6_25_2014ppb); 
meanIsoprene2015 = nanmean(subsetIsoprene2015ppb);
meanIsoprene2016 = nanmean(subsetIsoprene2016ppb); 


%%% Scatter Plot of 7/23/2009 Flight %%% 
figure
subplot(3,3,1)
scatter3(subsetLongitude7_23_2009, subsetLatitude7_23_2009, subsetAltitude7_23_2009,[80],subsetIsoprene7_23_2009,'filled')
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('7/23/2009 Isoprene ppt Levels','Fontsize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
axis([-119.9 -119.2 35.1 37.3])
hold off

%%% Scatter Plot of the 7/1/2010 %%%

subplot(3,3,2)
scatter3(subsetLongitude7_1_2010,subsetLatitude7_1_2010,subsetAltitude7_1_2010,[80],subsetIsoprene7_1_2010,'filled')
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('7/1/2010 Isoprene ppt Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])

%%% Scatter Plot of the 6/25/2014 Flight 
subplot(3,3,3)
scatter3(subsetLongitude6_25_2014,subsetLatitude6_25_2014,subsetAltitude6_25_2014,[80],subsetIsoprene6_25_2014,'filled')
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('6/25/2014 Isoprene ppt Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])
hold off 

%%% Scatter Plot 6/24/2015 Flight %%%
subplot(3,3,4)
scatter3(subsetLongitude2015,subsetLatitude2015,subsetAltitude2015,[80],subsetIsoprene2015,'filled');
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('6/24/2015 Isoprene ppt Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])
hold off 

%%% Scatter Plot 6/18/2016 Flight %%% 
subplot(3,3,5)
scatter3(subsetLongitude2016,subsetLatitude2016,subsetAltitude2016,[80],subsetIsoprene2016,'filled')
colormap jet
caxis([0,100])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('6/18/2016 Isoprene ppt Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])
hold off 

%%%% Begin the Vertical Distribution of Isoprene with Atltitude %%%%

figure
subplot(3,3,1);
scatter(subsetIsoprene7_23_2009ppb,subsetAltitude7_23_2009,[80],'filled'); 
xlabel('ppbv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Isoprene Vertical Distribution 7/23/2009','FontSize',18)

subplot(3,3,2);
scatter(subsetIsoprene7_1_2010ppb,subsetAltitude7_1_2010,[80],'filled')
xlabel('ppbv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Isoprene Vertical Distribution 7/1/2010','FontSize',18)

subplot(3,3,3);
scatter(subsetIsoprene6_25_2014ppb,subsetAltitude6_25_2014,[80],'filled')
xlabel('ppbv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Isoprene Vertical Distribution 6/25/2014','FontSize',18)

subplot(3,3,4);
scatter(subsetIsoprene2015ppb,subsetAltitude2015,[80],'filled')
xlabel('ppbv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Isoprene Vertical Distribution 6/24/2015','FontSize',18)

subplot(3,3,5);
scatter(subsetIsoprene2016ppb,subsetAltitude2016,[80],'filled')
xlabel('ppbv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Isoprene Vertical Distribution 6/18/2016','FontSize',18)




%%%% Begin Alpha Pinene Spatial Distribution Distribution %%%% 

% figure
% subplot(2,2,1)
% scatter3(subsetLongitude7_23_2009, subsetLatitude7_23_2009, subsetAltitude7_23_2009,[80],subsetAlphaPinene7_23_2009,'filled')
% colormap jet
% caxis([0,30])
% colorbar
% ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
% xlabel('Longitude', 'FontSize', 16)
% ylabel('Latitude','FontSize', 16)
% zlabel('Altitude (m)','FontSize', 16)
% title('7/23/2009 Alpha Pinene Levels','Fontsize',20)
% s=shaperead('CA_counties.shp','UseGeoCoords',true);
% hold on
% geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% axis([-119.9 -119.2 35.1 37.3])
% hold off

%%% AlphaPinene Scatter Plot of the 7/1/2010 %%%

subplot(2,2,1)
scatter3(subsetLongitude7_1_2010,subsetLatitude7_1_2010,subsetAltitude7_1_2010,[80],subsetAlphaPinene7_1_2010,'filled')
colormap jet
caxis([0,30])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('7/1/2010 Alpha Pinene Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])

% %%% Scatter Plot of the 6/25/2014 Flight 
% subplot(3,3,3)
% scatter3(subsetLongitude6_25_2014,subsetLatitude6_25_2014,subsetAltitude6_25_2014,[80],subsetIsoprene6_25_2014ppb,'filled')
% colormap jet
% caxis([0,1])
% colorbar
% ylabel(colorbar, 'ppbv','Interpreter', 'tex','FontSize',16)
% xlabel('Longitude', 'FontSize', 16)
% ylabel('Latitude','FontSize', 16)
% zlabel('Altitude (m)','FontSize', 16)
% title('6/25/2014 Isoprene ppb Levels','FontSize',20)
% s=shaperead('CA_counties.shp','UseGeoCoords',true);
% % s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
% hold on
% geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% % geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
% axis([-119.9 -119.2 35.1 37.3])
% hold off 

%%%  AlphaPinene Scatter Plot 6/24/2015 Flight %%%
subplot(2,2,2)
scatter3(subsetLongitude2015,subsetLatitude2015,subsetAltitude2015,[80],subsetAlphaPinene2015,'filled');
colormap jet
caxis([0,30])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('6/24/2015 Alpha Pinene Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])
hold off 
 
%%% AlphaPinene Scatter Plot 6/18/2016 Flight %%% 
subplot(2,2,3)
scatter3(subsetLongitude2016,subsetLatitude2016,subsetAltitude2016,[80],subsetAlphaPinene2016,'filled')
colormap jet
caxis([0,30])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex','FontSize',16)
xlabel('Longitude', 'FontSize', 16)
ylabel('Latitude','FontSize', 16)
zlabel('Altitude (m)','FontSize', 16)
title('6/18/2016 Alpha Pinene Levels','FontSize',20)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-119.9 -119.2 35.1 37.3])
hold off 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%% Begin Alpha Pinene Vertical Distribution %%%%% 

% figure
% subplot(2,2,1);
% scatter(subsetAlphaPinene7_23_2009,subsetAltitude7_23_2009,[80],'filled'); 
% xlabel('pptv','FontSize',16)
% ylabel('Altitude (m)','FontSize',16)
% title('Alpha Pinene Vertical Distribution 7/23/2009','FontSize',18)

subplot(2,2,1);
scatter(subsetAlphaPinene7_1_2010,subsetAltitude7_1_2010,[80],'filled')
xlabel('pptv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Alpha Pinene Vertical Distribution 7/1/2010','FontSize',18)

%%% Only one Alpha Pinene Data Point %%% 
% subplot(3,3,3);
% scatter(subsetAlphaPinene6_25_2014,subsetAltitude6_25_2014,[80],'filled')
% xlabel('pptv','FontSize',16)
% ylabel('Altitude (m)','FontSize',16)
% title('Alpha Pinene Vertical Distribution 6/25/2014','FontSize',18)

subplot(2,2,2);
scatter(subsetAlphaPinene2015,subsetAltitude2015,[80],'filled')
xlabel('pptv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Alpha Pinene Vertical Distribution 6/24/2015','FontSize',18)

subplot(2,2,3);
scatter(subsetAlphaPinene2016,subsetAltitude2016,[80],'filled')
xlabel('pptv','FontSize',16)
ylabel('Altitude (m)','FontSize',16)
title('Alpha Pinene Vertical Distribution 6/18/2016','FontSize',18)



%%% Max, Mean and Mins %%% 
Year = [2009 2010 2014 2015 2016];
MeanSubsetIsoprene7_23_2009 = nanmean(subsetIsoprene7_23_2009ppb); 
MeanSubsetIsoprene7_1_2010 = nanmean(subsetIsoprene7_1_2010ppb);
MeanSubsetIsoprene6_25_2014 = nanmean(subsetIsoprene6_25_2014ppb);
MeanSubsetIsoprene2015 = nanmean(subsetIsoprene2015ppb);
MeanSubsetIsoprene2016 = nanmean(subsetIsoprene2016ppb);

MeanIsopreneYearValues = [MeanSubsetIsoprene7_23_2009 MeanSubsetIsoprene7_1_2010 MeanSubsetIsoprene6_25_2014 MeanSubsetIsoprene2015 MeanSubsetIsoprene2016];
bar(Year,MeanIsopreneYearValues); 
ylabel('ppb','FontSize',25);
xlabel('Year','FontSize',25);
title('Mean Isoprene Levels for SARP Missions','FontSize',30); 

ts1 = timeseries(MeanIsopreneYearValues,Year)
plot(ts1)




%% Standard Deviations %%
StdSubsetIsoprene7_23_2009 = nanstd(subsetIsoprene7_23_2009ppb); 
StdubsetIsoprene7_1_2010 = nanstd(subsetIsoprene7_1_2010ppb);
StdSubsetIsoprene6_25_2014 = nanstd(subsetIsoprene6_25_2014ppb);
StdSubsetIsoprene2015 = nanstd(subsetIsoprene2015ppb);
StdSubsetIsoprene2016 = nanstd(subsetIsoprene2016ppb);

StdIsopreneYearValues = [nanstd(subsetIsoprene7_23_2009ppb) nanstd(subsetIsoprene7_1_2010ppb) nanstd(subsetIsoprene6_25_2014ppb) nanstd(subsetIsoprene2015ppb)  nanstd(subsetIsoprene2016ppb)];   





