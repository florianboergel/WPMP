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


x_vestas = (0:25)';
y_vestas = [0, 0,0,0,91,200,362,588,889,1255,1604,1769,1798,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800,1800]';

x_enercon = (0:25)';
y_enercon = [0,0,3,25,82,174,321,532,815,1180,1580,1900,2200,2400,2480,2700,2850,2950,3020,3020,3020,3020,3020,3020,3020,3020]';

for i = 1:25
    fino1_vestasyield(i,1) = y_vestas(i)*weibull_Fino1(i);
    fino2_vestasyield(i,1) = y_vestas(i)*weibull_Fino2(i);
    fino1_enerconyield(i,1) = y_enercon(i)*weibull_Fino1(i);
    fino2_enerconyield(i,1) = y_enercon(i)*weibull_Fino2(i);
end

if sum(fino1_vestasyield) > sum(fino2_vestasyield)
    disp('Fino 1 besser, vestas')
else
    disp('Fino 2 besser, vestas')
end

if sum(fino1_enerconyield) > sum(fino2_enerconyield)
    disp('Fino 1 besser, enercon')
else
    disp('Fino 2 besser, enercon')
end
 
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
windsOnHeight = [];
j = 1
for i = 1:length(fino1_v90)
    if (fino1_d90(i) >= 240 && fino1_d90(i) <= 285)
        wind_sector(j,2) = fino1_v90(1,i);
        wind_sector(j,1) = fino1_d90(1,i);
        windsOnHeight(j,1) = Fino1.ws33(1,i);
        windsOnHeight(j,2) = Fino1.ws40(1,i);
        windsOnHeight(j,3) = Fino1.ws50(1,i);
        windsOnHeight(j,4) = Fino1.ws60(1,i);
        windsOnHeight(j,5) = Fino1.ws70(1,i);
        windsOnHeight(j,6) = Fino1.ws80(1,i);
        windsOnHeight(j,7) = Fino1.ws90(1,i);
        windsOnHeight(j,8) = Fino1.ws100(1,i);
        j= j+1;
    end
end

wind_sector_mean = nanmean(wind_sector(:,2));
wind_prof = [];
for i = 1:100
    wind_prof(i,1) = 0.2/0.4 *(log(i/10^-6));
    wind_prof(i,2) = wind_sector_mean*(i/90)^(0.11);
end
figure();
hold on;
plot(1:100, wind_prof(:,1))
plot(1:100,wind_prof(:,2), 'o')
plot(90,wind_sector_mean,'*')
xlabel('Height in [m]');
ylabel('windspeed in [m/s]');
title('Vertical Profile');

avgPerHeight = [8];
for i=1:8
   avgPerHeight(i) = nanmean(windsOnHeight(:,i));
end
figure();
hold on;
plot([33,40,50,60,70,80,90,100],avgPerHeight(:), 'o')

logProfileModel = @(b,z) b(1)/b(2) *(log(z/b(3)));
empPowerModel = @(c,x) avgPerHeight(8)*((x/90).^c(1));
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
logProfileCoeffs = nlinfit([33,40,50,60,70,80,90,100],avgPerHeight,logProfileModel,[0.2,0.4,10^-6],opts);
fplot(@(z) logProfileCoeffs(1)/logProfileCoeffs(2) *(log(z/logProfileCoeffs(3))));
empPowerCoeff = nlinfit([33,40,50,60,70,80,90,100],avgPerHeight,empPowerModel,[0.11],opts);
fplot(@(z) avgPerHeight(8)*(z/90)^(empPowerCoeff));

xlabel('Height in [m]');
ylabel('windspeed in [m/s]');
title('Vertical Profile All Heights');
legend('Data', 'Log Profile','Empirical Power Profile');
saveas(gcf,'figures/verticalProfileFits.png')
hold off;

