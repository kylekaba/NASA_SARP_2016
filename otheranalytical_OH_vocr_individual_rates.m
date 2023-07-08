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
ylabel(colorbar, 'VOCR s^{-1}','Interpreter','Tex','FontSize',20)
xlabel('NO_x ppb','FontSize',20,'Interpreter','Tex')
ylabel('PO_3 ppb/hr','FontSize',20,'Interpreter','Tex')
title('Analytical O_3 Model',20,'Interpreter','Tex')

subplot(2,2,2)
scatter3(transectsubsetLongitudeNew,transectsubsetLatitudeNew,transectsubsetAltitudeNew,[40],vocr);
colormap jet
caxis([0 15])
colorbar
ylabel(colorbar,'VOCR s^{-1}','Interpreter','Tex','FontSize',20)
xlabel('Longitude,'FontSize',20,'Interpreter','Tex')
ylabel('Latitude','FontSize',20,'Interpreter','Tex')
zlabel('Altitude','FontSize',20,'Interpreter','Tex')
