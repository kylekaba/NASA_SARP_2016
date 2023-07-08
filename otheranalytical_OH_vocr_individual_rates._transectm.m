% analytical_OH written by JGM 031505
% need to make L HOx equal to P HOx for final balance
% now includes alkyl nitrate production - turn off by setting alpha = 0
% T=310;% hot AGU
% T=304.2;% mod
% T=298;% low
% tc = xx;
% T = xx+273;
T = transectsubsetTemperatureNew;
M= fun_air(transectsubsetTemperatureNew,transectsubsetPressureNew);
k1 = kohno2new(T,M);
k2 = kho2no(T,M); % effective ro2 or ho2 +no
k3 = kro2no(T,M);
ka = kho2ho2(T,M);
kb = kho2ro2(T,M);
kc = kro2ro2(T,M);

phox = transectsubsetpHOxtimestep; 
vocr = transectsubsetOHReacNew;
%n = 50;

noxppb= transectsubsetNOxppbnew;%[0.001:0.01:50];
nox=1E-9.*noxppb.*M; % coefficients below have been determined for NOx = 80% NO2, 20% NO

no=transectsubsetNOpptvnew.*0.001.*M.*(1E-9);;
no2=transectsubsetNO2pptvnew.*0.001.*M.*(1E-9);; % mod AGU

gamma=((vocr)./(k2.*no)).^2;
alpha=0.08; % branching ratio for AN formation

%quadratic formula solution for oh using equation PHOx=LHOx;
a = (2*ka + 2*kb + 2*kc).*(vocr./((1-alpha).*k2.*no)).^2;

b = k1.*no2 + (alpha.*k2.*vocr)./((1-alpha).*k2);

c = -phox;

oh=(-b+sqrt(b.^2-4.*a.*c))./(2.*a);

% using conservation of radicals PHOx = LHOx
ho2 = (vocr.*oh)./(k2.*no);
% final_ro2=2.*ho2; %(final ro2 is doubled to include ro2=ho2)

pozone=(k2+k3).*ho2.*no; %production rate in molec /cc s
pozone_ppbh=3600.*1E9.*pozone./M; %production rate in ppb/hr
phno3_ppbh=3600.*1E9.*k1.*oh.*no2./M; %production rate in ppb/hr
tau_ppbh=3600.*k1.*oh; %production rate in ppb/hr

% ppans_ppbh=3600.*1E9.*kpan_form(T,M).*(final_ro2./6).*no2./M; % production rate in ppb/hr
ohppt = oh.*1E12./M;
ho2ppt = ho2.*1E12./M;

hold on;
figure
subplot(2,2,1)
scatter(noxppb,pozone_ppbh,[40],vocr,'filled');
colormap jet
caxis([0,15])
colorbar
ylabel(colorbar, 'VOCR s^{-1}','Interpreter','Tex','FontSize',22)
xlabel('NO_x ppb','FontSize',22,'Interpreter','Tex')
ylabel('PO_3 ppb/hr','FontSize',22,'Interpreter','Tex')
title('Analytical O_3 Model: ','FontSize',28,'Interpreter','Tex')

subplot(2,2,2)
scatter3(transectsubsetLongitudeNew,transectsubsetLatitudeNew,transectsubsetAltitudeNew,[40],vocr);
colormap jet
caxis([0 15])
colorbar
ylabel(colorbar,'VOCR s^{-1}','Interpreter','Tex','FontSize',20)
xlabel('Longitude','FontSize',20,'Interpreter','Tex')
ylabel('Latitude','FontSize',20,'Interpreter','Tex')
zlabel('Altitude','FontSize',20,'Interpreter','Tex')
title('Analytical O_3 Model: Boundary Layer','Interpreter','Tex','FontSize',24)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
text(-119.01125833333333,35.44215555555555,'Bakersfield','FontSize',10);
text(-119.29212777777778,36.34723888888889, 'Visalia','FontSize',10);  
text(-119.77948333333333,36.767672222222224, 'Fresno','FontSize',10);
text(-119.64474166666668,36.34091388888889, 'Hanford','FontSize',10);
text(-120.47946388888889,37.31654444444444, 'Merced','FontSize',10);
text(-120.05923333333332,36.97851111111111, 'Madera','FontSize',10);
text(-118.24029722222222,34.06395555555555, 'Los Angeles','FontSize',10);
text(-119.17610833333333,34.21370277777778, 'Oxnard','FontSize',10);
text(-119.69840833333333,34.42899166666666, 'Santa Barbara','FontSize',10);
text(-120.99726111111111,37.66448611111111, 'Modesto','FontSize',10);
text(-121.29224166666667,37.982575000000004, 'Stockton','FontSize',10);
text(-120.77240833333333,38.352425000000004, 'Jackson','FontSize',10);
text(-120.67944722222222,38.20396388888889, 'San Andreas','FontSize',10);
text(-120.38189999999999,37.992380555555556, 'Sonora','FontSize',10);
text(-119.96443611111111,37.49040277777778, 'Mariposa','FontSize',10);
text(-117.28917777777778,34.11415277777778, 'San Bernadino','FontSize',10);
text(-117.39424166666667,33.96915555555556, 'Riverside','FontSize',10);
text(-117.79421944444444,33.710886111111115, 'Irvine','FontSize',10);
text(-117.15923611111111,32.74406666666667, 'San Diego','FontSize',10);
text(-118.19470277777778,33.770830555555555, 'Long Beach','FontSize',10);
text(-116.21716666666667,33.71935555555556, 'Indio','FontSize',10);
axis([-124 -114 33 39]);