%%% Looking at Various Biogenic Emissions from 6/18/2016 %%%

Latitude = KORUSAQHskpingdc820160618RA1(2:21539,3); 
Longitude = KORUSAQHskpingdc820160618RA1(2:21539,4); 
Altitude = KORUSAQHskpingdc820160618RA1(2:21539,8).*0.3048;

Isoprene(Isoprene == -99999) = NaN;
Isoprene(Isoprene < 0) = NaN; 
Isoprene2 = Isoprene;
Isoprene2ppb = Isoprene2.*0.001; 

% Creating the 3D Scatter Plot % 
scatter3(Longitude,Latitude,Altitude,[5],Isoprene2ppb);
colormap jet
caxis([0,.1])
colorbar
ylabel(colorbar, 'Isoprene (ppb)','Interpreter', 'tex')

title('6/18/2016 Isoprene Concentrations with Altitude', 'FontSize', 18, 'Interpreter', 'tex')
xlabel('Longitude', 'FontSize', 14)
ylabel('Latitude','FontSize', 14)
zlabel('Altitude (m)','FontSize', 14)


s=shaperead('CA_counties.shp','UseGeoCoords',true);
% s1=shaperead('MX_STATE.SHP', 'UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
% geoshow([s1.Lat],[s1.Lon],'Color','Black','LineWidth',2);
axis([-124 -114 33 42.5])


