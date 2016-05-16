%% Task 1
disp('Load Data')
load('WMP_WEnMet_data.mat');

fino1_v90 = Fino1.ws90;
fino1_d90 = Fino1.wd90;
fino2_v92 = Fino2.ws92;
fino2_d91 = Fino2.wd91;

% Plot
hold on;
figure();
wind_rose(fino1_d90,fino1_v90,'labtitle','Fino 1','percBg');
figure();
wind_rose(fino2_d91,fino2_v92,'labtitle','Fino 2','percBg');
hold off;

% Validate wind rose ...
direction = 175;
j=1;
speeds = NaN(length(fino1_v90),1);
for i = 1:length(fino1_v90)
    if (fino1_d90(i) >= direction && fino1_d90(i) <direction+10)
        speeds(j,1) = fino1_v90(1,i);
        j= j+1;
    end
end

count05 = speeds(speeds >0 & speeds <=5);
count510 = speeds(speeds >5 & speeds <=10);
count1015 = speeds(speeds >10 & speeds <=15);
count1520 = speeds(speeds >15 & speeds <=20);
count2025 = speeds(speeds >20 & speeds <=25);
count2530 = speeds(speeds >25 & speeds <=30);

freq05 = length(count05)/length(speeds(~isnan(speeds)));
freq510 = length(count510)/length(speeds(~isnan(speeds)));
freq1015 = length(count1015)/length(speeds(~isnan(speeds)));
freq1520 = length(count1520)/length(speeds(~isnan(speeds)));
freq2025 = length(count2025)/length(speeds(~isnan(speeds)));
%% Task 2
mean1 = nanmean(fino1_v90);
dev1 = nanstd(fino1_v90);

mean2 = nanmean(fino2_v92);
dev2 = nanstd(fino2_v92)

k_Fino1 = 1;
Func_Fino1 = @(k_Fino1) (mean1*mean1/(dev1*dev1))*((gamma(1+2/k_Fino1))/(gamma(1+1/k_Fino1))^2-1)-1 
k_Fino1 = fsolve(Func_Fino1,k_Fino1);
disp(k_Fino1);
A_Fino1 = mean1/gamma(1+1/k_Fino1)
weibull_Fino1 = wblpdf(1:30,A_Fino1,k_Fino1)

k_Fino2 = 1;
Func_Fino2 = @(k_Fino2) (mean2*mean2/(dev2*dev2))*((gamma(1+2/k_Fino2))/(gamma(1+1/k_Fino2))^2-1)-1 
k_Fino2 = fsolve(Func_Fino2,k_Fino2);
disp(k_Fino2);
A_Fino2 = mean2/gamma(1+1/k_Fino2)
weibull_Fino2 = wblpdf(1:30,A_Fino2,k_Fino2)


figure();
hold on;
histogram(fino1_v90, 'Normalization', 'pdf');
plot(weibull_Fino1)
xlabel('windspeed in [m/s]');
ylabel('Probability [%]');
title('Fino 1')
dim = [.7 .5 .3 .3];
annotation('textbox',dim,'String',{'k =',num2str(k_Fino1), 'A =' ,num2str(A_Fino1)},'FitBoxToText','on');
hold off;

figure();
hold on;
histogram(fino2_v92, 'Normalization', 'pdf');
plot(weibull_Fino2)
xlabel('windspeed in [m/s]');
ylabel('Probability [%]');
title('Fino 2')
dim = [.7 .5 .3 .3];
annotation('textbox',dim,'String',{'k =',num2str(k_Fino2), 'A =' ,num2str(A_Fino2)},'FitBoxToText','on');
hold off;



%% Task 3
wind_sector = [];
j = 1
for i = 1:length(fino1_v90)
    if (fino1_d90(i) >= 240 && fino1_d90(i) <= 285)
        wind_sector(j,2) = fino1_v90(1,i);
        wind_sector(j,1) = fino1_d90(1,i);
        j= j+1;
    end
end

wind_sector_mean = nanmean(wind_sector(:,2))



% disp('Loading Data ...')
% raw_data = readtable('1301.txt','Delimiter','tab');
% time_stamp = raw_data{:, {'Time'}};
% raw_data = raw_data{:, {'d90', 'd33', 'u100', 'u90', 'u80', ...
% 'u70', 'u60', 'u50', 'u40', 'u33'}};
% raw_data(raw_data==-999) = NaN;
% t = (datenum(time_stamp, 'yyyy-mm-dd HH:MM:SS'));
% clear time_stamp;
% t = round(t*24*3600);
% last_time = 31*24*3600;
% tnew=[1:1:last_time]';
% n = length(tnew);
% data_pp = NaN(length(tnew),10);
% for i = 1:length(raw_data(:,1))
%     data_pp(t(i)-t(1)+1, :) = raw_data(i, :);
% end
% time = (1:length(data_pp))';
% data_pp = [time, data_pp];
% save('data_pp.mat', 'data_pp', 'raw_data');
% clearvars;
% %% Evaluation
% load('data_pp.mat');
% 
% % means and 10 min means
% disp('Computing 10min means and stddev');  
% means_interval10 = NaN((length(data_pp(:,1))/600),20);
% for i = 1:length(data_pp(:,1))/600
%    %treat wind vanes
%    radians = data_pp((i-1)*600+1:i*600,2)/180*pi;
%    meanSin = nanmean(sin(radians));
%    meanCos = nanmean(cos(radians));
%    tanVal = atan2(meanSin,meanCos);
%    means_interval10(i,1) = tanVal*180/pi;
%    
%    radiansPrime = radians - tanVal;
%    primeUnwrapped = unwrap(radiansPrime);
%    unwrappedStddev = nanstd(primeUnwrapped);
%    means_interval10(i,2) = unwrappedStddev*180/pi;
% 
%    radians=data_pp((i-1)*600+1:i*600,3)/180*pi;
%    meanSin = nanmean(sin(radians));
%    meanCos = nanmean(cos(radians));
%    tanVal = atan2(meanSin,meanCos);
%    means_interval10(i,3) = tanVal*180/pi;
%    
%    radiansPrime = radians - tanVal;
%    primeUnwrapped = unwrap(radiansPrime);
%    unwrappedStddev = nanstd(primeUnwrapped);
%    means_interval10(i,4) = unwrappedStddev*180/pi;
%    %treat windspeeds u 
%    for j=3:10  
%        means_interval10(i,j*2-1) = nanmean(data_pp((i-1)*600+1:i*600,j+1));
%        means_interval10(i,j*2) = nanstd(data_pp((i-1)*600+1:i*600,j+1));
%    end
% end
% 
% save('meansAndStddev.mat', 'means_interval10');
% %% plot for January 30th
% load('meansAndStddev.mat');
% start30thJan = 29*24*6;
% figure();
% print -r1500
% plot(means_interval10(start30thJan+1:start30thJan+24*6,7), '-r');
% hold on;
% plot(means_interval10(start30thJan+1:start30thJan+24*6,7)+means_interval10(start30thJan+1:start30thJan+24*6,8), '--b');
% plot(means_interval10(start30thJan+1:start30thJan+24*6,7)-means_interval10(start30thJan+1:start30thJan+24*6,8), '--b');
% xlabel('10 minutes interval count [1]');
% ylabel('windspeed in [m/s]');
% legend('mean u_{90}','mean u_{90} \pm \sigma_u','Location','northwest');
% disp('saving plot to Plots/mean_interval_withstd.png');
% print('Plots/mean_interval_withstd.png','-dpng','-r500')
% hold off;
% clearvars;
% %% spikes
% load('meansAndStddev.mat');
% load('data_pp.mat');
% j= 1;
% k = 1;
% for i = 1:length(data_pp(:,1))/600
%     if max((data_pp((i-1)*600+1:i*600,4))) > (means_interval10(i,5)+5*means_interval10(i,6))
%         spikes(j,1) = i
%         j = j+1;
%     end
%     if min((data_pp((i-1)*600+1:i*600,4))) < (means_interval10(i,5)-5*means_interval10(i,6))
%         spikes(k,2) = i
%         k = k+1;
%     end
% end 
% 
% % Plot spikes
% i = 1084;
% window = 1 
% figure();
% for i = spikes(1:5,1)'
%     figure();
%     plot_data_mean(1:600,1) = means_interval10(i,5);
%     plot_data_mean(1:600,2) = means_interval10(i,5)+5*means_interval10(i,6);
%     plot_data_mean(1:600,3) = means_interval10(i,5)-5*means_interval10(i,6);
% 
%     hold on;
%     plot((data_pp((i-1)*600+1:i*600,4)), '-r')
%     plot(plot_data_mean(1:600,1),'-g')
%     plot(plot_data_mean(1:600,2),'-b')
%     plot(plot_data_mean(1:600,3),'-b')
%     xlabel('Time in [s]')
%     ylabel('Windspeed in [m/s]')
%     legend('Windspeed','Mean Windspeed','Standard Deviation', 'Location', 'northwest');
%     hold off;
%     stri = num2str(i);
%     str = strcat('Plots/spikesintervall', stri, '.png');
%     print(str,'-dpng','-r500');
% end
% disp('saving plot to Plots/10min_interval_with_highspikes.png')
% saveas(gcf,'Plots/10min_interval_with_highspikes.png')
% 
% window = 1 
% figure();
% for i = spikes(:,2)'
%     figure();
%     plot_data_mean(1:600,1) = means_interval10(i,5);
%     plot_data_mean(1:600,2) = means_interval10(i,5)+5*means_interval10(i,6);
%     plot_data_mean(1:600,3) = means_interval10(i,5)-5*means_interval10(i,6);
% 
%     hold on;
%     plot((data_pp((i-1)*600+1:i*600,4)), '-r')
%     plot(plot_data_mean(1:600,1),'-g')
%     plot(plot_data_mean(1:600,2),'-b')
%     plot(plot_data_mean(1:600,3),'-b')
%     xlabel('Time in [s]')
%     ylabel('Windspeed in [m/s]')
%     legend('Windspeed','Mean Windspeed','Standard Deviation', 'Location', 'northwest');
%     hold off;
%     stri = num2str(i);
%     str = strcat('Plots/spikesintervall', stri, '.png');
%     print(str,'-dpng','-r500');
% end
% %disp('saving plot to Plots/10min_interval_with_lowerspikes.png')
% %saveas(gcf,'Plots/10min_interval_with_lowerspikes.png')
% clearvars
% %% Task 6
% load('data_pp.mat');
% load('meansAndStddev.mat');
% tau = 1;
% for i = 1:length(data_pp(:,4))-1
%     du(i,1) = data_pp(i+tau,4)-data_pp(i,4);
% end
% 
% du = du/nanstd(du);
% [hist_y, hist_x] = hist(du,60);
% hist_y = hist_y/sum(hist_y);
% gausx = -10:1:10;
% gausy = (1/sqrt(2*pi))*exp(-gausx.^2/2);
% 
% figure();
% semilogy(hist_x, hist_y);
% hold on;
% semilogy(gausx,gausy);
% hold off;
% set(0,'DefaultTextInterpreter', 'latex');
% xlabel('$\delta u_\tau$');
% ylabel('$p(\delta u_\tau)$');
% print('Plots/tau_pdf_gauss10.png','-dpng','-r500')
% clearvars