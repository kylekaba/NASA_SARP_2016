s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',2);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
%axis([-110 -70 20 50])