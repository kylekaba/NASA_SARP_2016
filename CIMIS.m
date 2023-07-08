%%%% CIMIS Data %%%%

ETo = monthly(:,3); 
AvgSoilTemp = monthly(:,29);
MonthYear = monthly(:,4);

scatter(MonthYear, AvgSoilTemp)


startdate = datenum('01/2009')
enddate = datenum('06/2016')



x=linspace(startdate,enddate,41); 

startdate = datenum('01/2009')
enddate = datenum('06/2012')

%%% Organizing the Time %%%

datestr(73012,12)



%%% Xtick labels %%% 
plot(1:91,newEto,'LineWidth',2)
ax = gca;
ax.XTick = [7 19 30 42 54 66 78 90];
ax.XTickLabel = {monthly{7} monthly{19} monthly{30} monthly{42} monthly{54} monthly{66} monthly{78} monthly{90}}
ylabel('ETo','FontSize',19) 
title('ETo levels from Jan. 2009 - June 2016 Shafter Station') 


