% analytical_OH written by JGM 031505
% need to make L HOx equal to P HOx for final balance
% now includes alkyl nitrate production - turn off by setting alpha = 0
% T=310;% hot AGU
% T=304.2;% mod
% T=298;% low
% tc = xx;
% T = xx+273;

T= subsetTemperatureNew;%subsetTemperatureNew; %300;
M= fun_air(subsetTemperatureNew,subsetPressureNew); % air density % 
k1 = kohno2new(T,M);
k2 = kho2no(T,M); % effective ro2 or ho2 +no
k3 = kro2no(T,M);
ka = kho2ho2(T,M);
kb = kho2ro2(T,M);
kc = kro2ro2(T,M);


%% Total VOC Reactivity %% 
%% Athos OH - NOxReactivity 


%% Isoprene Reactivity with OH %%
% IsopreneOHReac = function_OrgReac(Isopreneppt,TemperatureNew,PressureNew,(2.7.*10^-11.*exp(390./TemperatureNew))); 
%% Note that NOxReac and NOx total in the workspace are the same %%% 

%subsetOHReacNew(subsetOHReacNew < 0) = NaN;

%OHandNOxDifference = subsetOHReacNew - subsetNOxReacnew; 

% Angelique's PHox 

subsetIsopreneOHReac(subsetIsopreneOHReac < 0) = NaN; 

phox = subsetpHOxtimestep; %2*0.5*M*1e-12; %% Setting my PHox to 0 to see the effects of drought % 
phox_pptv = (subsetpHOxtimestep./M).*1E12; 
vocr = subsetOHReacNew; %- IsopreneOHReac;  %2*20; %% Athos OH reactivity minus NOx reactivity %%

vocr(vocr < 0) = NaN;  
% 
% vocr_noisoprene = OHandNOxDifference - subsetIsopreneOHReac; 
% 
% vocr_noisopene(vocr_noisoprene < 0) = NaN; 


%%% Potential way to filter out the VOCrs %%%


%% Isoprene is 0.04% of the BVOC emissions %% 
% n = 50; % many times it cycles through % 

%noxppb = 0.5.*subsetNOxppbnew; %[0.001:0.01:50]; % Concentration of NO + NO2 ppb % 
%noxppb = 0.5.*subsetNOxppbnew; 

%% Converting to molecules/cm^3
nox = 1E-9.*noxppb.*M; % coefficients below have been determined for NOx = 80% NO2, 20% NO

% of the total NOx concentration, 25% is NO and 75% is NO2 %% Use real
% measured ratios, i.e. NO/NOx, NO2/NOx
no = 0.5.*subsetNOpptvnew.*0.001.*M.*(1E-9); %0.25.*nox; 

no2= 0.5*subsetNO2pptvnew.*0.001.*M.*(1E-9); %(0.75).*nox; % mod AGU

gamma=((vocr)./(k2.*no)).^2;
alpha = 0.02; % branching ratio for AN formation %Ange's has 0.02 for alpha

%quadratic formula solution for oh using equation PHOx=LHOx;
a = (2*ka + 2*kb + 2*kc).*(vocr./((1-alpha).*k2.*no)).^2;

b = k1.*no2 + (alpha.*k2.*vocr)./((1-alpha).*k2);

c = -phox;

oh=(-b+sqrt(b.^2-4.*a.*c))./(2.*a);

% using conservation of radicals PHOx = LHOx
ho2 = (vocr.*oh)./(k2.*no);
% final_ro2=2.*ho2; %(final ro2 is doubled to include ro2=ho2)

%% Making an Equation for Ozone Loss %% 
% Vd = Depsosition Velocity in cm/s %
%% Plugging in Test Values %% 
% Vd = 0.25;
% leafarea = 2;
% % LAI = 2;
% 
% 
% Volume_forppb = (1.66*10^-15*8.314.*Temperature*10^6)./Pressure; 
% O3_molecm3 = O3ppbv_2./Volume_forppb; 
% O3_molecm3(O3_molecm3 < 0) = NaN; 
% 
% lozone_dry = ((O3ppbv_2.*M.*Vd.*.3600.*leafarea));
% 


pozone =(k2+k3).*ho2.*no; %production rate in molec /cc s
pozone_ppbh= 3600.*1E9.*pozone./M; %production rate in ppb/hr 
% Above variable is with 50% NOx %
%pozone_ppbh_actualPO3 = 3600.*1E9.*pozone./M;

pozone_ppbh_Isoprene60per = 3600.*1E9.*pozone./M;
phno3_ppbh=3600.*1E9.*k1.*oh.*no2./M; %production rate in ppb/hr
tau_ppbh=3600.*k1.*oh; %production rate in ppb/hr

%% Difference of 100% Nox - 50% Nox multiplied by 6 (6 hour production)%
pozonedifference_8hour = (pozone_ppbh_actualPO3 - pozone_ppbh).*8;

% difference of 60 % Isoprene% 
pozonedifference_Isoprene60per = (pozone_ppbh_actualPO3 - pozone_ppbh_Isoprene60per).*8; 
% ppans_ppbh=3600.*1E9.*kpan_form(T,M).*(final_ro2./6).*no2./M; % production rate in ppb/hr
ohppt = oh.*1E12./M;
ho2ppt = ho2.*1E12./M;

hold on;
subplot(2,2,1)
figure
scatter(noxppb,pozone_ppbh,[40],subsetTemperatureNew,'filled') %ohppt% , 'm--')
colormap jet
caxis([290,300])
colorbar
ylabel(colorbar, 'Temperature','Interpreter','Tex','FontSize',20)
xlabel('NO_x ppb','FontSize',20,'Interpreter','Tex')
ylabel('PO_3 ppb/hr','FontSize',20,'Interpreter','Tex')
title('Analytical O_3 Model: Boundary Layer','Interpreter','Tex','FontSize',24); 
ylim([0 25] )


%% VOCR Scatter Plot %% 
subplot(2,2,2);
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],subsetTemperatureNew,'filled')
colormap jet
caxis([290,300])
colorbar
ylabel(colorbar, 'Temperature','Interpreter','Tex','FontSize',18)
title('6/18/2016 Temperature Boundary Layer Spatial Plot', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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




hold on;
%subplot(2,2,3)
figure 
scatter(noxppb,pozone_ppbh,[40],vocr,'filled') %ohppt% , 'm--')
colormap jet
caxis([0,10])
colorbar
ylabel(colorbar, 'VOCR (s^{-1})','Interpreter','Tex','FontSize',20)
xlabel('NO_x ppb','FontSize',20,'Interpreter','Tex')
ylabel('PO_3 ppb/hr','FontSize',20,'Interpreter','Tex')
title('Analytical O_3 Model: SJV Boundary Layer','Interpreter','Tex','FontSize',24); 
ylim([0 25] )


%% VOCR 3D Scatter Plot %% 
% subplot(2,2,4);
figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],vocr,'filled')
colormap jet
caxis([0 15])
colorbar
ylabel(colorbar, 'VOCR (s^{-1})','Interpreter','Tex','FontSize',18)
title('6/18/2016 VOCR Boundary Layer Spatial Plot', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
text(-119.01125833333333,35.44215555555555,'Bakersfield','FontSize',10);
text(-119.29212777777778,36.34723888888889, 'Visalia','FontSize',10);  
text(-119.77948333333333,36.767672222222224, 'Fresno','FontSize',10);
text(-119.64474166666668,36.34091388888889, 'Hanford','FontSize',10);
text(-120.47946388888889,37.31654444444444, 'Merced','FontSize',10);
text(-120.05923333333332,36.97851111111111, 'Madera','FontSize',10);
%text(-118.24029722222222,34.06395555555555, 'Los Angeles','FontSize',10);
%text(-119.17610833333333,34.21370277777778, 'Oxnard','FontSize',10);
%text(-119.69840833333333,34.42899166666666, 'Santa Barbara','FontSize',10);
%text(-120.99726111111111,37.66448611111111, 'Modesto','FontSize',10);
%text(-121.29224166666667,37.982575000000004, 'Stockton','FontSize',10);
%text(-120.77240833333333,38.352425000000004, 'Jackson','FontSize',10);
%text(-120.67944722222222,38.20396388888889, 'San Andreas','FontSize',10);
%text(-120.38189999999999,37.992380555555556, 'Sonora','FontSize',10);
%text(-119.96443611111111,37.49040277777778, 'Mariposa','FontSize',10);
%text(-117.28917777777778,34.11415277777778, 'San Bernadino','FontSize',10);
%text(-117.39424166666667,33.96915555555556, 'Riverside','FontSize',10);
%text(-117.79421944444444,33.710886111111115, 'Irvine','FontSize',10);
%text(-117.15923611111111,32.74406666666667, 'San Diego','FontSize',10);
%text(-118.19470277777778,33.770830555555555, 'Long Beach','FontSize',10);
%text(-116.21716666666667,33.71935555555556, 'Indio','FontSize',10);
axis([-124 -114 33 39]);


%%% Time to for loop This Mother %%% 

tick = 1;
for i = 1:length(AltitudeNew); % all points in 2009 scatterplot
    if AltitudeNew(i) < 1200 && LatitudeNew(i) > 35 && LatitudeNew(i) < 37.5 && LongitudeNew(i) < -118.5 && LongitudeNew(i) > -121
%         subsetLatitudeNew(tick) = LatitudeNew(i);
%         subsetLongitudeNew(tick) = LongitudeNew(i);
%         subsetIsopreneOHReac(tick) = IsopreneOHReac(i);
%         subsetAltitudeNew(tick) = AltitudeNew(i); 
%         subsetNO2pptvnew(tick) = NO2pptvnew(i); 
%         subsetNOpptvnew(tick) = NOpptvnew(i);
%         subsetOHReacNew(tick) = OHReacNew(i);
%         subsetpHOxtimestep(tick) = pHOxtimestep(i);
%         subsetTemperatureNew(tick) = TemperatureNew(i);
%         subsetPressureNew(tick) = PressureNew(i);
%         subsetNOxReacnew(tick) = NOxtotal(i);
%         subsetIsopreneOHReac(tick) = IsopreneOHReac(i); 
         subsetNOxppbnew(tick) = NOxppbnew(i); 
        % subsetO3ppbv_2(tick) = O3ppbv_2(i);
        subsetUTC(tick) = UTC(i);
        
        tick = tick + 1;
    end
end 
       
figure
plot(subsetTemperatureNew,vocr);

%%%% Ozone Color Map %%%% 
figure
%subplot(2,2,1)
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],pozone_ppbh,'filled')
colormap jet
caxis([0 15])
colorbar
ylabel(colorbar, 'PO_3 ppb/hr','Interpreter','Tex','FontSize',18)
title('6/18/2016 Modeled PO_3 Production', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
text(-119.01125833333333,35.44215555555555,'Bakersfield','FontSize',11);
text(-119.29212777777778,36.34723888888889, 'Visalia','FontSize',11);  
text(-119.77948333333333,36.767672222222224, 'Fresno','FontSize',11);
text(-119.64474166666668,36.34091388888889, 'Hanford','FontSize',11);
text(-120.47946388888889,37.31654444444444, 'Merced','FontSize',11);
text(-120.05923333333332,36.97851111111111, 'Madera','FontSize',11);
text(-118.24029722222222,34.06395555555555, 'Los Angeles','FontSize',11);
text(-119.17610833333333,34.21370277777778, 'Oxnard','FontSize',11);
text(-119.69840833333333,34.42899166666666, 'Santa Barbara','FontSize',11);
text(-120.99726111111111,37.66448611111111, 'Modesto','FontSize',11);
text(-121.29224166666667,37.982575000000004, 'Stockton','FontSize',11);
text(-120.77240833333333,38.352425000000004, 'Jackson','FontSize',11);
text(-120.67944722222222,38.20396388888889, 'San Andreas','FontSize',11);
text(-120.38189999999999,37.992380555555556, 'Sonora','FontSize',11);
text(-119.96443611111111,37.49040277777778, 'Mariposa','FontSize',11);
text(-117.28917777777778,34.11415277777778, 'San Bernadino','FontSize',11);
text(-117.39424166666667,33.96915555555556, 'Riverside','FontSize',11);
text(-117.79421944444444,33.710886111111115, 'Irvine','FontSize',11);
text(-117.15923611111111,32.74406666666667, 'San Diego','FontSize',11);
text(-118.19470277777778,33.770830555555555, 'Long Beach','FontSize',11);
text(-116.21716666666667,33.71935555555556, 'Indio','FontSize',11);
axis([-124 -114 33 39]);

%% Section on color coding the plot with pHOx %% 

% Analytical Model of PO3 
figure 
scatter(noxppb,pozone_ppbh,[40],phox_pptv,'filled') %ohppt% , 'm--')
colormap jet
caxis([0,0.6])
colorbar
ylabel(colorbar, 'pHO_x (pptv)','Interpreter','Tex','FontSize',20)
xlabel('NO_x ppb','FontSize',20,'Interpreter','Tex')
ylabel('PO_3 ppb/hr','FontSize',20,'Interpreter','Tex')
title('Analytical O_3 Model: SJV Boundary Layer','Interpreter','Tex','FontSize',24); 
ylim([0 25] )





%PHOx Spatial MATLAB Plot 
figure
%subplot(2,2,1)
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],phox_pptv,'filled')
colormap jet
caxis([0 0.5])
colorbar
ylabel(colorbar, 'pHO_x','Interpreter','Tex','FontSize',18)
title('6/18/2016 Modeled PO_3 Production', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
s=shaperead('CA_counties.shp','UseGeoCoords',true);
hold on
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',1);
text(-119.01125833333333,35.44215555555555,'Bakersfield','FontSize',11);
text(-119.29212777777778,36.34723888888889, 'Visalia','FontSize',11);  
text(-119.77948333333333,36.767672222222224, 'Fresno','FontSize',11);
text(-119.64474166666668,36.34091388888889, 'Hanford','FontSize',11);
text(-120.47946388888889,37.31654444444444, 'Merced','FontSize',11);
text(-120.05923333333332,36.97851111111111, 'Madera','FontSize',11);
text(-118.24029722222222,34.06395555555555, 'Los Angeles','FontSize',11);
text(-119.17610833333333,34.21370277777778, 'Oxnard','FontSize',11);
text(-119.69840833333333,34.42899166666666, 'Santa Barbara','FontSize',11);
text(-120.99726111111111,37.66448611111111, 'Modesto','FontSize',11);
text(-121.29224166666667,37.982575000000004, 'Stockton','FontSize',11);
text(-120.77240833333333,38.352425000000004, 'Jackson','FontSize',11);
text(-120.67944722222222,38.20396388888889, 'San Andreas','FontSize',11);
text(-120.38189999999999,37.992380555555556, 'Sonora','FontSize',11);
text(-119.96443611111111,37.49040277777778, 'Mariposa','FontSize',11);
text(-117.28917777777778,34.11415277777778, 'San Bernadino','FontSize',11);
text(-117.39424166666667,33.96915555555556, 'Riverside','FontSize',11);
text(-117.79421944444444,33.710886111111115, 'Irvine','FontSize',11);
text(-117.15923611111111,32.74406666666667, 'San Diego','FontSize',11);
text(-118.19470277777778,33.770830555555555, 'Long Beach','FontSize',11);
text(-116.21716666666667,33.71935555555556, 'Indio','FontSize',11);
axis([-124 -114 33 39]);

%% 
%%% NOx Color Plot %%% 
figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],subsetNOxppbnew,'filled')
colormap jet
caxis([0 10])
colorbar
ylabel(colorbar, 'NO_x ppb','Interpreter','Tex','FontSize',18)
title('6/18/2016 NO_x Color Map ', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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


%% Test Control Halving NOx by 0.5. 

figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],subsetNOxppbnew,'filled')
colormap jet
caxis([0 10])
colorbar
ylabel(colorbar, 'NO_x ppb','Interpreter','Tex','FontSize',18)
title('6/18/2016 NO_x Color Map ', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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

%% Actual O3 Production in the SJV, not calculated %%
%figure
subplot(2,2,2)
figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],noxppb_2,'filled')
colormap jet
caxis([0 3])
colorbar
ylabel(colorbar, 'NO_x ppb','Interpreter','Tex','FontSize',18)
title('6/18/2016 NO_x Color Map ', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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


%% Plotting the Scatter plot of measured O3 levels %%
subplot(2,2,2)
figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],subsetO3ppbv_2new,'filled')
colormap jet
caxis([30 80])
colorbar
ylabel(colorbar, 'O_3 ppb','Interpreter','Tex','FontSize',18)
title('6/18/2016 Measured O_3 Color Map ', 'FontSize', 28)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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


%%% Trying to Find the Tranasect Region %%% 
figure
scatter3(subsetLongitudeNew,subsetLatitudeNew,subsetAltitudeNew,[40],subsetUTC); 
colormap jet
caxis([62000 80000])
colorbar
ylabel(colorbar, 'O_3 ppb','Interpreter','Tex','FontSize',18)
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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




%% 
%%%% Extracting Transect Regions %%%% 

tick = 1;
for i = 1:length(AltitudeNew); % all points in 2009 scatterplot
    if AltitudeNew(i) < 1200 && LatitudeNew(i) > 35 && LatitudeNew(i) < 37.5 && LongitudeNew(i) < -118.5 && LongitudeNew(i) > -121 && UTC(i) < 69000 && UTC(i) > 67000
        transectsubsetLatitudeNew(tick) = LatitudeNew(i);
        transectsubsetLongitudeNew(tick) = LongitudeNew(i);
        transectsubsetIsopreneOHReac(tick) = IsopreneOHReac(i);
        transectsubsetAltitudeNew(tick) = AltitudeNew(i); 
        transectsubsetNO2pptvnew(tick) = NO2pptvnew(i); 
        transectsubsetNOpptvnew(tick) = NOpptvnew(i);
        transectsubsetOHReacNew(tick) = OHReacNew(i);
        transectsubsetpHOxtimestep(tick) = pHOxtimestep(i);
        transectsubsetTemperatureNew(tick) = TemperatureNew(i);
        transectsubsetPressureNew(tick) = PressureNew(i);
        transectsubsetNOxReacnew(tick) = NOxtotal(i);
        transectsubsetIsopreneOHReac(tick) = IsopreneOHReac(i); 
        transectsubsetNOxppbnew(tick) = NOxppbnew(i); 
        transectsubsetO3ppbv_2(tick) = O3ppbv_2(i);
        transectsubsetUTC(tick) = UTC(i);
        
       
        tick = tick + 1;
    end
end 

%% Storing the New Variables underneath %%

transectsubsetLatitudeNew = transpose(transectsubsetLatitudeNew);
transectsubsetLongitudeNew = transpose(transectsubsetLongitudeNew);  
transectsubsetIsopreneOHReac = transpose(transectsubsetIsopreneOHReac); 
transectsubsetAltitudeNew = transpose(transectsubsetAltitudeNew);
transectsubsetNO2pptvnew = transpose(transectsubsetNO2pptvnew);
transectsubsetNOpptvnew = transpose(transectsubsetNOpptvnew);
transectsubsetOHReacNew = transpose(transectsubsetOHReacNew);
transectsubsetpHOxtimestep = transpose(transectsubsetpHOxtimestep);
transectsubsetTemperatureNew = transpose(transectsubsetTemperatureNew);
transectsubsetPressureNew = transpose(transectsubsetPressureNew);
transectsubsetNOxReacnew = transpose(transectsubsetNOxReacnew);
transectsubsetIsopreneOHReac = transpose(transectsubsetIsopreneOHReac);
transectsubsetNOxppbnew = transpose(transectsubsetNOxppbnew);
transectsubsetO3ppbv_2 = transpose(transectsubsetO3ppbv_2);
transectsubsetUTC = transpose( transectsubsetUTC); 




%%

figure
scatter3(transpose(transectsubsetLongitudeNew),transpose(transectsubsetLatitudeNew),transpose(transectsubsetAltitudeNew)); 
xlabel('Longitude', 'FontSize', 18)
ylabel('Latitude','FontSize', 18)
zlabel('Altitude (m)','FontSize', 18)
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
