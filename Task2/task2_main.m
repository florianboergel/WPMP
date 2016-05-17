%% Task 1
disp('Load Data')
load('WMP_WEnMet_data.mat');

fino1_v90 = Fino1.ws90;
fino1_d90 = Fino1.wd90;
fino2_v92 = Fino2.ws92;
fino2_d91 = Fino2.wd91;

% Plot
WindRose(fino1_d90,fino1_v90,'AngleNorth',0,'AngleEast',90);
saveas(gcf,'figures/WindRose_Fino1.png')
WindRose(fino2_d91,fino2_v92,'AngleNorth',0,'AngleEast',90);
saveas(gcf,'figures/WindRose_Fino2.png')


% Validate wind rose ...
[a,b] = ksdensity(fino2_d91)
figure();
plot(b,a)
title('Wind Direction density fino2')
saveas(gcf,'figures/Validation_WindRose_Fino2.png')
[c,d] = ksdensity(fino1_d90)
figure();
plot(d,c)
title('Wind Direction density fino1')
saveas(gcf,'figures/Validation_WindRose_Fino1.png')
%% Task 2
mean1 = nanmean(fino1_v90);
dev1 = nanstd(fino1_v90);

mean2 = nanmean(fino2_v92);
dev2 = nanstd(fino2_v92);

% interpolate
k_Fino1 = 1;
Func_Fino1 = @(k_Fino1) (mean1*mean1/(dev1*dev1))*((gamma(1+2/k_Fino1))/(gamma(1+1/k_Fino1))^2-1)-1 
k_Fino1 = fsolve(Func_Fino1,k_Fino1);
disp(k_Fino1);
A_Fino1 = mean1/gamma(1+1/k_Fino1);
weibull_Fino1 = wblpdf(1:30,A_Fino1,k_Fino1);

k_Fino2 = 1;
Func_Fino2 = @(k_Fino2) (mean2*mean2/(dev2*dev2))*((gamma(1+2/k_Fino2))/(gamma(1+1/k_Fino2))^2-1)-1 
k_Fino2 = fsolve(Func_Fino2,k_Fino2);
disp(k_Fino2);
A_Fino2 = mean2/gamma(1+1/k_Fino2);
weibull_Fino2 = wblpdf(1:30,A_Fino2,k_Fino2);


figure();
hold on;
histogram(fino1_v90, 'Normalization', 'pdf');
plot(weibull_Fino1)
xlabel('windspeed in [m/s]');
ylabel('Probability [%]');
title('Fino 1')
dim = [.7 .5 .3 .3];
annotation('textbox',dim,'String',{'k =',num2str(k_Fino1), 'A =' ,num2str(A_Fino1)},'FitBoxToText','on');
saveas(gcf,'figures/Hist_withfit_Fino1.png')
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
saveas(gcf,'figures/Hist_withfit_Fino2.png')
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