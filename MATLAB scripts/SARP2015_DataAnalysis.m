%% Looking at the SARP 2015 Merge File %%

Latitude2015 = SARP2015June24Flight(:,7);
Longitude2015 = SARP2015June24Flight(:,8);
Altitude2015 = SARP2015June24Flight(:,12);


%% Isoprene in pptv %% 
Isoprene2015 = SARP2015June24Flight(:,74);
Isoprene2015(Isoprene2015 == -888) = NaN;
Isoprene2015(Isoprene2015 < 0) = NaN; 
Isoprene2015_2 = Isoprene2015; 

scatter3(Longitude2015,Latitude2015,Altitude2015,[40],Isoprene2015_2,'filled');
colormap jet
caxis([0,1000])
colorbar
ylabel(colorbar, 'pptv','Interpreter', 'tex')
xlabel('Longitude', 'FontSize', 14)
ylabel('Latitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)
title('SARP 6/24/2015 Flight Isoprene Levels','FontSize',24) 
s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 42.5])

%%%% Now Looking at SARP 2014 Data %%%%

%% SARP Flight on 6/24/2014 %% 
Isoprene6_24_14 = June242014Flight1(:,28);
Isoprene6_24_14(Isoprene6_24_14 == -888) = NaN;
Isoprene6_24_14(Isoprene6_24_14 < 0) = NaN;
Isoprene6_24_2014_2 = Isoprene6_24_14; 


%%% 6/25/2015 %%%
Isoprene6_25_14 = June252014(:,28);
Isoprene6_25_14(Isoprene6_25_14 == -888) = NaN;
Isoprene6_25_14(Isoprene6_25_14 < 0) = NaN;
IsopreneIsoprene6_25_14_2 = Isoprene6_25_14; 


