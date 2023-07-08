%%% figure
pcolor(Longitude,Latitude,July2010);
shading interp
hold on
contour(Longitude,Latitude,July2010);
text(-119.5,36.746841,'Fresno','FontSize',12)
text(-118.5,35.5,'Bakersfield','FontSize',12)
geoshow(36.746841,-119.772591,'DisplayType', 'Point','MarkerEdgeColor','k','Linewidth',2,'MarkerSize',11)
geoshow(35.3733,-119.0187,'DisplayType', 'Point','MarkerEdgeColor','k','Linewidth',2,'MarkerSize',11)
geoshow([s.Lat],[s.Lon],'Color','Black','LineWidth',2);
geoshow('SARP_02.shp','Color', 'r', 'MarkerEdgeColor', 'auto')
ylim([34 39])
xlim([-124 -112])
title('Isoprene Emissions')
xlabel('Longitude')
ylabel('Latitude')
colorbar
hold off 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% for loop %%%%
%% Trying to extract data from this stupid file %%

A = zeros(721,1441);
A(1,2:1441) = Longitude; 
A(2:721,1) = Latitude; 
A(2:721,2:1441) = July2010; 


%%% Getting LAI Data From July 2009 %%% 
tick_lat = 1;
tick_lon = 1;
LAI_2009 = zeros(40,40);
for i = 1:720 ; 
   for  j = 1:1440;
    if Longitude(j) >= -124 && Longitude(j) <= -114 && Latitude(i) <= 42 && Latitude(i) >= 32
        LAI_2009(tick_lat,tick_lon) = July2009(i,j);  
        tick_lon = tick_lon + 1;
        if tick_lon == 41 ;
            tick_lon = 1 ;
            tick_lat = tick_lat + 1;
        end 
        end 
    end
end

LAI_2009(LAI_2009 == 99999) = NaN; 

sum(sum(LAI_2009); 


%%% Attempting to Determine The Amount of NaNs in the DataSet %%% 
isnan(LAI_2009);

NANs2009 = sum(sum(isnan(LAI_2009)))
        
        