%% Loading

close all; clc;
load('201401020600_WLS200S-17_data_.mat');

ts_VAD = data.ts;
az_c_VAD = data.az_c;
az_c_VAD = az_c_VAD + 180
rs_VAD = data.rs;
el_VAD = data.el;
rg_VAD = data.rg;
x_VAD = data.x;
y_VAD = data.y;
z_VAD = data.z;
cnr_VAD = data.cnr;
%% Step 1 
% Write Data do one double struct
data_VAD = [ts_VAD az_c_VAD rs_VAD cnr_VAD el_VAD rg_VAD x_VAD y_VAD z_VAD];

% Filter data by CnR, all DATA?
cnr_ind = find(cnr_VAD <= -20 | cnr_VAD >= -5);
%data_VAD(cnr_ind, 2:3) = NaN;
%rs_VAD(cnr_ind) = NaN;
%az_VAD(cnr_ind) = NaN;


% get all range gates
count = 1

for q = 1:100
    a = find(rg_VAD == q);
    if length(a) ~= 0 
        rg_ind(:,count) = a;
        count = count + 1;
    end    
end

% sort by range gates from lowest to highest
for i = 1:14
    rs_VAD_filter(:,i) = data_VAD(rg_ind(:,i),3);
    az_c_VAD_filter(:,i) = data_VAD(rg_ind(:,i),2);
    ts_VAD_filter(:,i) = data_VAD(rg_ind(:,i),1);
end

% how long does one full 360° scan take?
time_per_scan = 360 / 25
time_per_scan_intervalls = time_per_scan / 0.4

%% Step 2
phi = az_c_VAD_filter(500:1000,1);
v_r = rs_VAD_filter(500:1000,1);

VADCos = @(param,phi) param(1)*cosd(phi-param(2))+param(3)
startvalues = [0, 0, 0]
fitparam = lsqcurvefit(VADCos, startvalues, phi, v_r)


