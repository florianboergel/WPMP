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
startvalues = [0, 0, 0];
numberOfScans = floor(length(az_c_VAD_filter)/intervalls_per_scan-1)
calc_times = zeros(numberOfScans,14);
calc_vHor = zeros(numberOfScans,14);
calc_vVer = zeros(numberOfScans,14);
calc_D = zeros(numberOfScans,14);

for j = 1:14
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
            v_ver = -fitparam(3) / cosd(60);
            D=fitparam(1)+180;
            calc_times(i,j) = ts_VAD_filter((i-1)*intervalls_per_scan+18,j);
            calc_vHor(i,j) =  v_hor;
            calc_vVer(i,j) = v_ver;
            calc_D(i,j) = D;
        end
    end 
end


%% Step 3: compare lidar and sonic data
%10 minute averages of the lidar data 
startTime = datenum('02-Jan-2014 06:00:00','dd-mmm-yyyy HH:MM:SS');

for j=1:14
    for interval = 1:36 %we have 36 ten minute intervals
        avgs_horSpeed(interval,j) = mean(calc_vHor(calc_times(:,j) >= startTime +(interval-1)/(24*6) & ... 
                                                    calc_times(:,j) < startTime +interval/(24*6),j));
        avgs_horDir(interval,j) = mean(calc_D(calc_times(:,j) >= startTime +(interval-1)/(24*6) & ... 
                                                calc_times(:,j) < startTime +interval/(24*6),j));
    end;
end;

%calc vertical profiles and compare



