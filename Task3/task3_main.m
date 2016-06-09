%% Task 1
load('WMP_WEnMet_data.mat');
raw_data = readtable('MERRA2_N55.000_E013.125_1992-2016.txt','Delimiter','tab','HeaderLines',26);

fino2_v92 = Fino2.ws92';
fino2_d91 = Fino2.wd91';

%compute indices for relevant timestamps in both formats (Merra2 and Fino2)
timestamps = raw_data.Var1(:,1);
first_timestamp = find(strcmp(timestamps(:), '01.01.2010 00:00'));
last_timestamp = find(strcmp(timestamps(:), '31.05.2014 23:00'));
last_timestampFino2 = find(Fino2.time==datenum('31-May-2014 23:55:00'));

% Check status flags of Merra2 Data
if any(raw_data.Var9(:)) ~= 0
    disp('error');
end

%creating 1h averages for FINO2 over relevant period
fino2_1h_92 = [];
for i = 1:last_timestampFino2/6
    fino2_1h_92(i,1) = datenum('01-Jan-2010 00:00:00')+(i-1)*1/24;
    range_array_v = fino2_v92((i-1)*6+1:i*6,1);
    range_array_dir = fino2_d91((i-1)*6+1:i*6,1);
    fino2_1h_92(i,3) = nanmean(range_array_v);
    fino2_1h_92(i,2) = nanmean(range_array_dir);
end

%connected_data contains 1h averages of merra2 and fino2 data as well as
%corresponding directions and timestamps
connected_data(:,1) = fino2_1h_92(:,1); %timstamps (1h intervals from 01.01.2010 until 31.05.2014)
connected_data(:,2) = fino2_1h_92(:,3); %speed fino 2
connected_data(:,3) = fino2_1h_92(:,2); %direction fino2
connected_data(:,4) = raw_data.Var2(first_timestamp:last_timestamp); %speed merra2
connected_data(:,5) = raw_data.Var3(first_timestamp:last_timestamp); %direction merra2

for i = 1:length(connected_data(:,1))
    if raw_data.Var8(i) ~= 0 || raw_data.Var9(i) ~= 0 || isnan(connected_data(i,2)) == 1 || isnan(connected_data(i,3)) == 1
        connected_data(i,2:5) = NaN;
    end
end


%% Task 2
for i=1:39
    sortedCell{i} = [];
end;

%compute sector index and put the relevant data (merra2 and fino2 speed and
%timestamp) into cell structure by appending
for i = 1:length(connected_data(:,1))
    if ~isnan(connected_data(i,3)) && ~isnan(connected_data(i,4)) && ~isnan(connected_data(i,2))
        sortIndex = floor(connected_data(i,3)/360*12)+1;
        sortedCell{sortIndex*3-2} = [sortedCell{sortIndex*3-2}, connected_data(i,1)]; %timestamp
        sortedCell{sortIndex*3-1} = [sortedCell{sortIndex*3-1}, connected_data(i,4)]; %merra 2
        sortedCell{sortIndex*3} = [sortedCell{sortIndex*3}, connected_data(i,2)]; %fino 2
        sortedCell{37} = [sortedCell{37}, connected_data(i,1)]; %timestamp
        sortedCell{38} = [sortedCell{38}, connected_data(i,4)]; %merra 2
        sortedCell{39} = [sortedCell{39}, connected_data(i,2)]; %fino 2
    end;
end;

%% Task 3
avgSectorPerMonth = zeros(53,27);
firstMonth = datenum('01-Jan-2010');

%first column: month indices
for i=1:53
    avgSectorPerMonth(i,1) = addtodate(firstMonth,i-1,'month');
end;

%second column: monthly averages sector independent
for i=1:53
     valuesInMonthAndSectorMerra2 = sortedCell{38}(find((sortedCell{37} > addtodate(firstMonth,i-1,'month')) ...
            & (sortedCell{37} <= addtodate(firstMonth,i,'month'))));
     valuesInMonthAndSectorFino2 = sortedCell{39}(find((sortedCell{37} > addtodate(firstMonth,i-1,'month')) ...
            & (sortedCell{37} <= addtodate(firstMonth,i,'month'))));
     avgSectorPerMonth(i,2)= nanmean(valuesInMonthAndSectorMerra2);
     avgSectorPerMonth(i,3)= nanmean(valuesInMonthAndSectorFino2);
end;

%%next: wind speeds for both locations as monthly average per sector
for sectorIndex = 1:12
    for monthIndex = 1:53
        valuesInMonthAndSectorMerra2 = sortedCell{sectorIndex*3-1}(find((sortedCell{sectorIndex*3-2} > addtodate(firstMonth,monthIndex-1,'month')) ...
            & (sortedCell{sectorIndex*3-2} <= addtodate(firstMonth,monthIndex,'month'))));
       valuesInMonthAndSectorFino2 = sortedCell{sectorIndex*3}(find((sortedCell{sectorIndex*3-2} > addtodate(firstMonth,monthIndex-1,'month')) ...
            & (sortedCell{sectorIndex*3-2} <= addtodate(firstMonth,monthIndex,'month'))));
        avgSectorPerMonth(monthIndex,sectorIndex*2+2)= nanmean(valuesInMonthAndSectorMerra2);
        avgSectorPerMonth(monthIndex,sectorIndex*2+3)= nanmean(valuesInMonthAndSectorFino2);
    end;
end;

%% Task 4 Task 5
regressionParameters = zeros(13,4);

figure();
hold on;
[regressionParameters(1,2),  regressionParameters(1,3), regressionParameters(1,4)] = ... 
        regression(avgSectorPerMonth(:,2), avgSectorPerMonth(:,3), 'one');
regressionParameters(1,1) = regressionParameters(1,2) * regressionParameters(1,2);
fplot( @(x) regressionParameters(1,3)*x + regressionParameters(1,4), [0 18]);
scatter(avgSectorPerMonth(:,2), avgSectorPerMonth(:,3));
xlabel('MERRA-2 speed [m/s]','FontSize',12);
ylabel('FINO2 speed [m/s]','FontSize',12);
title('All Sectors', 'FontSize',20);
legend(strcat('R^2= ',num2str(regressionParameters(1,1)),'  m= ',num2str(regressionParameters(1,3)),' b= ',num2str(regressionParameters(1,4))),'Data points','Location','northwest');
saveas(gcf,strcat('figures/scatterPlot_AllSectors.jpg'));
hold off;

for i=2:13
    figure();
    hold on;
    [regressionParameters(i,2),  regressionParameters(i,3), regressionParameters(i,4)] = ... 
        regression(avgSectorPerMonth(:,i*2), avgSectorPerMonth(:,i*2+1), 'one');
    regressionParameters(i,1) = regressionParameters(i,2) * regressionParameters(i,2);
    fplot( @(x) regressionParameters(i,3)*x + regressionParameters(i,4),[0 max(avgSectorPerMonth(:,i*2))]);
    scatter(avgSectorPerMonth(:,i*2), avgSectorPerMonth(:,i*2+1));
    xlabel('MERRA-2 speed [m/s]','FontSize',12);
    ylabel('FINO2 speed [m/s]','FontSize',12);
    title(strcat(num2str((i-2)*30), '° - ', num2str((i-1)*30) ,'°'),'FontSize',20);
    legend(strcat('R^2= ',num2str(regressionParameters(i,1)),'  m= ',num2str(regressionParameters(i,3)),' b= ',num2str(regressionParameters(i,4))),'Data points','Location','northwest');
    saveas(gcf,strcat('figures/scatterPlot_Sector',num2str((i-2)*30),'.jpg'));
    hold off;
end;


%% Task 6
for i = 1:length(raw_data.Var2(:))
    sectorIndex= floor(raw_data.Var3(i) /360*12)+1;
    longTermCorrectedFino2(i,1) = raw_data.Var2(i) * regressionParameters(sectorIndex,2) + regressionParameters(sectorIndex,3);
    longTermCorrectedFino2(i,2) = raw_data.Var3(i);
end;

%% Task 7
weibullParam = zeros(13,2);
weibullParam(1,:) = wblfit(longTermCorrectedFino2(:,1));

for sectorIndex=1:12
    weibullParam(sectorIndex+1,:) = wblfit(longTermCorrectedFino2(find(longTermCorrectedFino2(:,2) >= (sectorIndex-1)*30 & ... 
        longTermCorrectedFino2(:,2) < sectorIndex*30),1));
end;

%% Task 8
figure();
hold on;
plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam(1,1), weibullParam(1,2)),'--','LineWidth',2);
for sectorIndex=2:13
    plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam(sectorIndex,1), weibullParam(sectorIndex,2)),'LineWidth',2);
end;
xlabel('Windspeed in [m/s]');
ylabel('Prob.');
title('Weibull distribution sector-wise');
legend('All', 'Sector 1', 'Sector 2', 'Sector 3', 'Sector 4', 'Sector 5', 'Sector 6', 'Sector 7', ...
        'Sector 8', 'Sector 9', 'Sector 10', 'Sector 11', 'Sector 12','Location','northeast');
saveas(gcf,'figures/weibullPlots.jpg');
hold off; 

%% Task 9
weibullParam_shortTerm = zeros(13,2);
weibullParam_shortTerm(1,:) = wblfit(sortedCell{39});
figure('visible','off');
hold on;
plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam(1,1), weibullParam(1,2)),'LineWidth',2);
plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam_shortTerm(1,1),weibullParam_shortTerm(1,2)),'LineWidth',2);
xlabel('Windspeed in [m/s]');
ylabel('Prob.');
title('All Sectors','FontSize',20);
legend('Long-Term Corrected','Short-Term Measured', 'Location','northeast');
saveas(gcf,'figures/shortVsLongTerm_AllSectors.jpg');
hold off; 
    
for sectorIndex=2:13
    weibullParam_shortTerm(sectorIndex,:) = wblfit(sortedCell{(sectorIndex-1)*3});
    figure();
    hold on;
    plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam(sectorIndex,1), weibullParam(sectorIndex,2)),'LineWidth',2);
    plot(0:0.1:25,wblpdf(0:0.1:25,weibullParam_shortTerm(sectorIndex,1),weibullParam_shortTerm(sectorIndex,2)),'LineWidth',2);
    xlabel('Windspeed in [m/s]');
    ylabel('Prob.');
    title(strcat(num2str((sectorIndex-2)*30), '° - ', num2str((sectorIndex-1)*30) ,'°'),'FontSize',20);
    legend('Long-Term Corrected','Short-Term Measured', 'Location','northeast');
    saveas(gcf,strcat('figures/shortVsLongTerm_Sector',num2str((sectorIndex-2)*30),'.jpg'));
    hold off; 
end;

%% Task 10
overviewTable = zeros(3,13);
for i=1:13
    overviewTable(1,i) = regressionParameters(i,1);
    overviewTable(2,i) = weibullParam(i,1) / weibullParam_shortTerm(i,1) ;
    overviewTable(3,i) = weibullParam(i,2) / weibullParam_shortTerm(i,2) ;
end;

