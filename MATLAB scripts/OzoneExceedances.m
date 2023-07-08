%%%% Creating Ozone Exceedance Data %%%%

Year = [2010 2011 2012 2013 2014 2015];
Fresno = [40 56 48 34 70 56];
Merced = [31 41 25 31 44 34];
Bakersfield = [48 51 83 47 39 54]; 
Visalia = [57 33 60 10 27 52];

figure
bar(Year,Fresno);
xlabel('Year', 'FontSize',16)
ylabel('Days Exceeding EPA 8-hr O_3 standard','Interpreter','Tex','FontSize',16)
title('Fresno','FontSize',18)

figure
bar(Year,Merced);
xlabel('Year', 'FontSize',16)
ylabel('Days Exceeding EPA 8-hr O_3 standard','Interpreter','Tex','FontSize',16)
title('Merced','FontSize',18)

figure
bar(Year,Bakersfield); 
xlabel('Year', 'FontSize',16)
ylabel('Days Exceeding EPA 8-hr O_3 standard','Interpreter','Tex','FontSize',16)
title('Bakersfield','FontSize',18)

figure
bar(Year,Visalia); 
xlabel('Year', 'FontSize',16)
ylabel('Days Exceeding EPA 8-hr O_3 standard','Interpreter','Tex','FontSize',16)
title('Visalia','FontSize',18)