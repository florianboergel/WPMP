%% Task 1: Import data and preprocess
disp('Loading Data ...')
raw_data = readtable('1301.txt','Delimiter','tab');
time_stamp = raw_data{:, {'Time'}};
raw_data = raw_data{:, {'d90', 'd33', 'u100', 'u90', 'u80', ...
'u70', 'u60', 'u50', 'u40', 'u33'}};
disp('Converting time with datenum ...')
t = (datenum(time_stamp, 'yyyy-mm-dd HH:MM:SS'));
clear time_stamp;
t = round(t*24*3600);

%% Task 2: Marking invalid data
disp('Marking invalid data ...')
%raw_data(raw_data==-999) = NaN(size(raw_data(raw_data==-999)));
raw_data(raw_data==-999) = NaN;

%% Task 3: Generate continuous time axis
disp('Creating continous time axis')
tnew=[t(1):1:t(end)]';
n = length(tnew);
data_pp = NaN(length(tnew),10);
disp('Writing preprocessed Data...')
for i = 1:length(raw_data(:,1))
    data_pp(t(i)-t(1)+1, :) = raw_data(i, :);
end
time = (1:length(data_pp))';
data_pp = [time, data_pp];
    save('data_pp.mat', 'data_pp', 'raw_data');

%% Task 4: 10 minutes mean, stddev and plotting
disp('Computing 10min means and stddev');
load('data_pp.mat');
means_interval10 = zeros(round(n/600),20);
for i = 1:n/600
   %treat wind vanes
   radians = data_pp((i-1)*600+1:i*600,2)/180*pi;
   meanSin = nanmean(sin(radians));
   meanCos = nanmean(cos(radians));
   tanVal = atan2(meanSin,meanCos);
   means_interval10(i,1) = tanVal*180/pi;
   
   radiansPrime = radians - tanVal;
   primeUnwrapped = unwrap(radiansPrime);
   unwrappedStddev = nanstd(primeUnwrapped);
   means_interval10(i,2) = unwrappedStddev*180/pi;

   radians=data_pp((i-1)*600+1:i*600,3)/180*pi;
   meanSin = nanmean(sin(radians));
   meanCos = nanmean(cos(radians));
   tanVal = atan2(meanSin,meanCos);
   means_interval10(i,3) = tanVal*180/pi;
   
   radiansPrime = radians - tanVal;
   primeUnwrapped = unwrap(radiansPrime);
   unwrappedStddev = nanstd(primeUnwrapped);
   means_interval10(i,4) = unwrappedStddev*180/pi;
   %treat windspeeds u 
   for j=3:10  
       means_interval10(i,j*2-1) = nanmean(data_pp((i-1)*600+1:i*600,j+1));
       means_interval10(i,j*2) = nanstd(data_pp((i-1)*600+1:i*600,j+1));
   end
end
save('meansAndStddev.mat', 'means_interval10');

%plot for January 30th u90
start30thJan = 29*24*6;
figure();
plot(means_interval10(start30thJan+1:start30thJan+24*6,7), '-r');
hold on;
plot(means_interval10(start30thJan+1:start30thJan+24*6,7)+means_interval10(start30thJan+1:start30thJan+24*6,8), '--b');
plot(means_interval10(start30thJan+1:start30thJan+24*6,7)-means_interval10(start30thJan+1:start30thJan+24*6,8), '--b');
xlabel('interval number');
ylabel('windspeed in [m/s]');
legend('Mean windspeed','Standard Deviation');
disp('saving plot to Plots/mean_intervall_withstd.png');
saveas(gcf,'Plots/mean_intervall_withstd.png');
hold off;

%% Task 5: Finding spike for u100
found = 0;
for i = 1:n/600
   currentMean = means_interval10(i,5);
   currentStddev = means_interval10(i,6);
   for j = 1:600
        if abs(data_pp((i-1)*600+j,4) - currentMean) > currentStddev*5
            %found spike in this interval
            disp(['found spike in interval ', num2str(i), ' with mean ', num2str(currentMean), ' and stddev ', num2str(currentStddev)]);
            interval = i;
            found = 1;
            break;
        end
   end
   if found == 1
       break;
   end
end
%plot the found interval containing a spike
figure();
plot(data_pp((i-1)*600+1:i*600,4), '-r');
meanCurve = @(x) currentMean;
upper = @(x) currentMean + currentStddev;
lower = @(x) currentMean - currentStddev;
upper5 = @(x) currentMean + 5*currentStddev;
lower5 = @(x) currentMean - 5*currentStddev;
hold on;
fplot(meanCurve,[0 600],'k-o');
fplot(upper,[0 600],'g--');
fplot(lower,[0 600],'g--');
fplot(upper5,[0 600],'r-.');
fplot(lower5,[0 600],'r-.');
xlabel('time in [s]');
ylabel('windspeed in [m/s]');
legend('Data','Mean windspeed','Standard Deviation','5 x Standard Deviation');
disp('saving plot to Plots/interval_spike.png');
saveas(gcf,'Plots/interval_spike.png');
hold off;

%% Task 6: Increment PDF

