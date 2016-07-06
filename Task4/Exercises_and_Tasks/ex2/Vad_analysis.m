%% Loading

close all; clc;
load('201401020600_WLS200S-17_data_.mat');

ts_VAD = data.ts;
az_c_VAD = data.az_c;
az_c_VAD = az_c_VAD + 180;
rs_VAD = data.rs;
el_VAD = data.el;
rg_VAD = data.rg;
x_VAD = data.x;
y_VAD = data.y;
z_VAD = data.z;
cnr_VAD = data.cnr;

%% Step 1 
% Write Data to one double struct
data_VAD = [ts_VAD az_c_VAD rs_VAD cnr_VAD el_VAD rg_VAD x_VAD y_VAD z_VAD];

% Filter data by CnR, by setting to NaN and later in Step 2 discarding NaN
% values
cnr_ind = find(cnr_VAD <= -20 | cnr_VAD >= -5);
rs_VAD(cnr_ind) = NaN;
az_c_VAD(cnr_ind) = NaN;


% get all range gates
count = 1;

for q = 1:100
    a = find(rg_VAD == q);
    if length(a) ~= 0 
        rg_ind(:,count) = a;
        rangeGateSet(count)=q;
        count = count + 1;
    end    
end

% separate data by range gates
for i = 1:14
    rs_VAD_filter(:,i) = data_VAD(rg_ind(:,i),3);
    az_c_VAD_filter(:,i) = data_VAD(rg_ind(:,i),2);
    ts_VAD_filter(:,i) = data_VAD(rg_ind(:,i),1);
end

% how long does one full 360° scan take?
time_per_scan = 360 / 25;
intervalls_per_scan = time_per_scan / 0.4;

%% Step 2
VADCos = @(param,phi) param(1)*cosd(phi-param(2))+param(3);
startvalues = [2, 10, 2];
numberOfScans = floor(length(az_c_VAD_filter)/intervalls_per_scan-1)
calc_times = zeros(numberOfScans,14);
calc_vHor = zeros(numberOfScans,14);
calc_vVer = zeros(numberOfScans,14);
calc_D = zeros(numberOfScans,14);

for j = 1:14 % number of Range Gates
    for i = 1:numberOfScans
        phi = az_c_VAD_filter((i-1)*intervalls_per_scan+18:(i-1)*intervalls_per_scan+53,j);
        v_r = rs_VAD_filter((i-1)*intervalls_per_scan+18:(i-1)*intervalls_per_scan+53,j);
        nanIndices = isnan(phi) | isnan(v_r); 
        phi(nanIndices) = [];
        v_r(nanIndices) = [];
        if length(phi>0)
            fitparam = lsqcurvefit(VADCos, startvalues, phi, v_r);
            savedparam{i,j} = fitparam;
            %calculate horiz and vertical speed
            v_hor = fitparam(1) / cosd(60);
            v_ver = -fitparam(3) / sind(60);
            D=(fitparam(1)+180);
            calc_times(i,j) = ts_VAD_filter((i-1)*intervalls_per_scan+18,j);
            calc_vHor(i,j) =  v_hor;
            calc_vVer(i,j) = v_ver;
            calc_D(i,j) = D;
        end
    end 
end


%% Step 3: compare lidar and sonic data
startTime = datenum('02-Jan-2014 06:00:00','dd-mmm-yyyy HH:MM:SS');
endTime = datenum('02-Jan-2014 12:00:00','dd-mmm-yyyy HH:MM:SS');

%10 minute averages of the lidar wind speed (horizontal) 
for j=1:14
    for interval = 1:36 %we have 36 ten minute intervals
        tenMinAvg_horSpeed(interval,j) = nanmean(calc_vHor(calc_times(:,j) >= startTime +(interval-1)/(24*6) & ... 
                                                    calc_times(:,j) < startTime +interval/(24*6),j));
        tenMinAvg_direction(interval,j) = nanmean(calc_D(calc_times(:,j) >= startTime +(interval-1)/(24*6) & ... 
                                                    calc_times(:,j) < startTime +interval/(24*6),j));
    end;
end;

%load fino1 data
fino_33 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_33m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_40 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_40m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_50 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_50m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_60 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_60m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_70 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_70m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_80 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_80m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_90 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_90m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_100 = readtable('fino-150504160459/FINO1_Windgeschwindigkeit_100m_20131220_20140121.dat','Delimiter','tab','HeaderLines',8);
fino_speeds(:,1) = fino_33.Var2((datenum(fino_33.Var1(:)) >= startTime & datenum(fino_33.Var1(:)) < endTime));
fino_speeds(:,2) = fino_40.Var2((datenum(fino_40.Var1(:)) >= startTime & datenum(fino_40.Var1(:)) < endTime));
fino_speeds(:,3) = fino_50.Var2((datenum(fino_50.Var1(:)) >= startTime & datenum(fino_50.Var1(:)) < endTime));
fino_speeds(:,4) = fino_60.Var2((datenum(fino_60.Var1(:)) >= startTime & datenum(fino_60.Var1(:)) < endTime));
fino_speeds(:,5) = fino_70.Var2((datenum(fino_70.Var1(:)) >= startTime & datenum(fino_70.Var1(:)) < endTime));
fino_speeds(:,6) = fino_80.Var2((datenum(fino_80.Var1(:)) >= startTime & datenum(fino_80.Var1(:)) < endTime));
fino_speeds(:,7) = fino_90.Var2((datenum(fino_90.Var1(:)) >= startTime & datenum(fino_90.Var1(:)) < endTime));
fino_speeds(:,8) = fino_100.Var2((datenum(fino_100.Var1(:)) >= startTime & datenum(fino_100.Var1(:)) < endTime));
fino_speeds(fino_speeds==-999) = NaN;
%overall mean speeds for lidar
for j=1:14
    avg_horSpeed_Lidar(j) = nanmean(tenMinAvg_horSpeed(:,j));
end;
%overall mean speeds for fino
avg_horSpeed_Fino(1) = nanmean(fino_speeds(:,1));
avg_horSpeed_Fino(2) = nanmean(fino_speeds(:,2));
avg_horSpeed_Fino(3) = nanmean(fino_speeds(:,3));
avg_horSpeed_Fino(4) = nanmean(fino_speeds(:,4));
avg_horSpeed_Fino(5) = nanmean(fino_speeds(:,5));
avg_horSpeed_Fino(6) = nanmean(fino_speeds(:,6));
avg_horSpeed_Fino(7) = nanmean(fino_speeds(:,7));
avg_horSpeed_Fino(8) = nanmean(fino_speeds(:,8));


%perform power log fit for fino1
logProfileModel = @(b,z) b(1)/0.4 *(log(z/b(2)));
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
logProfileCoeffsFino = real(nlinfit([33,40,50,60,70,80,90,100],avg_horSpeed_Fino,logProfileModel,[0.1,10^-5],opts));
lidarHeights = rangeGateSet*sind(60);
logProfileCoeffsLidar= real(nlinfit(lidarHeights,avg_horSpeed_Lidar,logProfileModel,[0.1,10^-5],opts));
%%
%plot both vertical wind speed profiles
x_range = [33,40,50,60,70,80,90,100]
figure();
hold on;
[xLogFino,yLogFino]=fplot(@(z) logProfileCoeffsFino(1)/0.4 *(log(z/logProfileCoeffsFino(2))),[0 100]);
[xLogLidar,yLogLidar]=fplot(@(z) logProfileCoeffsLidar(1)/0.4 *(log(z/logProfileCoeffsLidar(2))),[0 100]);
plot(yLogFino,xLogFino,'Color','b');
plot(yLogLidar,xLogLidar,'Color','r');
plot(avg_horSpeed_Lidar,lidarHeights,'or');
plot(avg_horSpeed_Fino,x_range, 'ob');
ylabel('Height in [m]');
xlabel('windspeed in [m/s]');
title('Vertical Profiles');
legend('Fino1','LiDAR VAD','Location','northwest');
saveas(gcf,'figures/verticalProfiles.png')
hold off;

%% Winddirections
fino_90_dir = readtable('fino-150504160459/FINO1_Windrichtung_90m_20131220_20140121.dat','Delimiter','tab','HeaderLines',6);
fino_90_dir_interval = fino_90_dir.Var2((datenum(fino_90_dir.Var1(:)) >= startTime & datenum(fino_90_dir.Var1(:)) < endTime));

WindRose(fino_90_dir_interval,fino_speeds(:,7),'AngleNorth',0,'AngleEast',90,'freqlabelangle',45,'MaxFrequency',6);
saveas(gcf,'figures/WindRose_Fino1.png')
WindRose(tenMinAvg_direction(:,7), tenMinAvg_horSpeed(:,7),'AngleNorth',0,'AngleEast',90,'freqlabelangle',45,'MaxFrequency',6);
saveas(gcf,'figures/WindRose_lidar.png')
figure();
hold on;
plot(fino_90_dir_interval);
plot(tenMinAvg_direction);
hold off;


%% Plots for report
figure();
hold on;
plot(VADCos(savedparam{1,1},1:360));
scatter(az_c_VAD_filter(18:53,1), rs_VAD_filter(18:53,1));
hold off;
xlabel('azimuth angle [degree]')
ylabel('V_{los} [m/s]')
saveas(gcf,'figures/fit.png')