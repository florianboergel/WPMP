%% Task 1
disp('Load Data')
load('WMP_WEnMet_data.mat');
raw_data = readtable('MERRA2_N55.000_E013.125_1992-2016.txt','Delimiter','tab','HeaderLines',26);

fino2_v92 = Fino2.ws92';
fino2_d91 = Fino2.wd91';

%compute timestamps
timestamps = raw_data.Var1(:,1);
first_timestamp = find(strcmp(timestamps(:), '01.01.2010 00:00'));
last_timestamp = find(strcmp(timestamps(:), '31.05.2014 23:00'));
last_timestampFino2 = find(Fino2.time==datenum('31-May-2014 23:55:00'));

% Check for continous time axis
count = 0
for i = 2:length(fino2_v92)
    if Fino2.time(i)-Fino2.time(i-1) >= 1/24/6*1.5
        count = i;
    end
end

if any(raw_data.Var9(:)) ~= 0
    disp('error');
end

fino2_1h_v92 = [];
for i = 1:last_timestampFino2/6
    fino2_1h_v92(i,1) = datenum('01-Jan-2010 00:00:00')+(i-1)*1/24;
    range_array_v = fino2_v92((i-1)*6+1:i*6,1);
    range_array_dir = fino2_d91((i-1)*6+1:i*6,1);

    fino2_1h_v92(i,3) = nanmean(range_array_v);
    fino2_1h_v92(i,2) = nanmean(range_array_dir);
end

connected_data(:,1) = fino2_1h_v92(:,1);
connected_data(:,2) = fino2_1h_v92(:,3); %speed fino 2
connected_data(:,3) = fino2_1h_v92(:,2); %direction fino2
connected_data(:,4) = raw_data.Var2(first_timestamp:last_timestamp); %speed merra2
connected_data(:,5) = raw_data.Var3(first_timestamp:last_timestamp); %direction merra2

count = 0 
for i = 1:length(connected_data(:,1))
    if raw_data.Var8(i) ~= 0 || raw_data.Var9(i) ~= 0 || isnan(connected_data(i,2)) == 1 || isnan(connected_data(i,3)) == 1
        connected_data(i,2:5) = NaN;
        count = count +1;
    end
end


%% Task 2
for i=1:36
    sortedCell{i} = [];
end;

for i = 1:length(connected_data(:,1))
    if ~isnan(connected_data(i,3))
        sortIndex = floor(connected_data(i,3)/360*12)+1;
        sortedCell{sortIndex*3-2} = [sortedCell{sortIndex*3-2}, connected_data(i,1)]; %timestamp
        sortedCell{sortIndex*3-1} = [sortedCell{sortIndex*3-1}, connected_data(i,4)]; %merra 2
        sortedCell{sortIndex*3} = [sortedCell{sortIndex*3}, connected_data(i,2)]; %fino 2
    end;
end;

%% Task 3
avgSectorPerMonth = zeros(53,25);
firstMonth = datenum('01-Jan-2010');

for i=1:53
    avgSectorPerMonth(i,1) = addtodate(firstMonth,i-1,'month');
end;

for sectorIndex = 1:12
    for monthIndex = 1:53
        valuesInMonthAndSectorMerra2 = sortedCell{sectorIndex*3-1}(find((sortedCell{sectorIndex*3-2} > addtodate(firstMonth,monthIndex-1,'month')) ...
            & (sortedCell{sectorIndex*3-2} <= addtodate(firstMonth,monthIndex,'month'))));
       valuesInMonthAndSectorFino2 = sortedCell{sectorIndex*3}(find((sortedCell{sectorIndex*3-2} > addtodate(firstMonth,monthIndex-1,'month')) ...
            & (sortedCell{sectorIndex*3-2} <= addtodate(firstMonth,monthIndex,'month'))));
        avgSectorPerMonth(monthIndex,sectorIndex*2)= nanmean(valuesInMonthAndSectorMerra2);
        avgSectorPerMonth(monthIndex,sectorIndex*2+1)= nanmean(valuesInMonthAndSectorFino2);
    end;
end;

%% Task 4 Task 5
regressionParameters = zeros(12,4);
for i=1:12
    subplot(3,4,i);
    %scatter(avgSectorPerMonth(:,i*2), avgSectorPerMonth(:,i*2+1));
    [regressionParameters(i,2),  regressionParameters(i,3), regressionParameters(i,4)] =  regression(avgSectorPerMonth(:,i*2), avgSectorPerMonth(:,i*2+1), 'one');
    regressionParameters(i,1) = regressionParameters(i,2)*regressionParameters(i,2);

    plotregression(avgSectorPerMonth(:,i*2), avgSectorPerMonth(:,i*2+1), 'Regression');
   % title(strcat('Sector  ', num2str((i-1)*30), '° - ', num2str(i*30) ,'°'));
    saveas(gcf,strcat('figures/scatterPlot_Sector',num2str((i-1)*30),'.jpg'));

end;


%% Task 6
for i = 1:length(raw_data.Var2(:))
    sectorIndex= floor(raw_data.Var3(i) /360*12)+1;
    %%slow: longTermCorrectedFino2(i,1) = datenum(raw_data.Var1(i));
    longTermCorrectedFino2(i,1) = raw_data.Var2(i) *regressionParameters(sectorIndex,2) + regressionParameters(sectorIndex,3);
    longTermCorrectedFino2(i,2) = raw_data.Var3(i);
end;

%% Task 7
weibullParam = zeros(13,2);
weibullParam(1,:) = wblfit(longTermCorrectedFino2(:,1));

for sectorIndex=1:12
    weibullParam(sectorIndex+1,:) = wblfit(longTermCorrectedFino2(find(longTermCorrectedFino2(:,2) >= (sectorIndex-1)*30 & ... 
        longTermCorrectedFino2(:,2) < sectorIndex*30),1));
%     figure;
%     hold on;
%     plot(wblpdf(1:30,weibullParam(sectorIndex+1,1), weibullParam(sectorIndex+1,2)));
%     histogram(longTermCorrectedFino2(find(longTermCorrectedFino2(:,2) >= (sectorIndex-1)*30 & ... 
%         longTermCorrectedFino2(:,2) < sectorIndex*30),1), 'Normalization', 'pdf');
%     hold off; 
end;

%% Task 8
figure();
hold on;
    plot(wblpdf(1:0.1:25,weibullParam(1,1), weibullParam(1,2)),'--','LineWidth',2);
for sectorIndex=2:13
    plot(wblpdf(1:0.1:25,weibullParam(sectorIndex,1), weibullParam(sectorIndex,2)),'LineWidth',2);
end;
xlabel('Windspeed in [m/s]');
ylabel('Prob.');
title('Weibull distribution sector-wise');
legend('All', 'Sector 1', 'Sector 2', 'Sector 3', 'Sector 4', 'Sector 5', 'Sector 6', 'Sector 7', ...
        'Sector 8', 'Sector 9', 'Sector 10', 'Sector 11', 'Sector 12','Location','northeast');
    saveas(gcf,'figures/weibullPlots.jpg');
hold off; 

%% Task 9
for sectorIndex=2:13
    shortTermParams = wblfit(sortedCell{(sectorIndex-1)*3});
    figure();
    hold on;
    plot(wblpdf(1:0.1:25,weibullParam(sectorIndex,1), weibullParam(sectorIndex,2)),'LineWidth',2);
    plot(wblpdf(1:0.1:25,shortTermParams(1),shortTermParams(2)),'LineWidth',2);
    xlabel('Windspeed in [m/s]');
    ylabel('Prob.');
    title('Short-term and long-term corrected data for Fino2');
    legend('Short-Term Measured', 'Long-Term Corrected','Location','northeast');
    saveas(gcf,strcat('figures/shortVsLongTerm_Sector',num2str((sectorIndex-2)*30),'.jpg'));
    hold off; 
end;

%% Task 10

