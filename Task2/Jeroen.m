%%
clear all;
close all;
clc;

%% Set constants
fs_LOS = 2.0;
fs_VAD = 2.5;

%% Load data
disp('Load data');
load('Data/FINO1.mat');
LOS_files = dir('Data/LOS');
VAD_files = dir('Data/VAD');

ts_VAD = [];
rs_VAD = [];
cnr_VAD = [];
az_VAD = [];

for i = 1:5 %length(VAD_files)
    disp(strcat('Loading VAD file ', num2str(i), ' of ', num2str(length(VAD_files))));
    if (VAD_files(i).isdir==0)
        load(strcat('Data/VAD/', VAD_files(i).name));
        ts_VAD = [ts_VAD; data.ts];
        rs_VAD = [rs_VAD; data.rs];
        cnr_VAD = [cnr_VAD; data.cnr];
        az_VAD = [az_VAD; data.az_c];
    end
end

az_VAD = az_VAD + 180;
data_VAD = [ts_VAD rs_VAD cnr_VAD az_VAD];

ts_LOS = [];
rs_LOS = [];
cnr_LOS = [];


for i = 1:5 %length(LOS_files)
    disp(strcat('Loading LOS file ', num2str(i), ' of ', num2str(length(LOS_files))));
    if (LOS_files(i).isdir==0)
        load(strcat('Data/LOS/', LOS_files(i).name));
        ts_LOS = [ts_LOS; data.ts];
        rs_LOS = [rs_LOS; data.rs];
        cnr_LOS = [cnr_LOS; data.cnr];
    end
end

data_LOS = [ts_LOS rs_LOS cnr_LOS];

%% Format time
disp('Format time');
starttime_VAD = min(data_VAD(:, 1));
data_VAD(:, 1) = data_VAD(:, 1) - min(data_VAD(:, 1));
data_VAD(:, 1) = data_VAD(:, 1) * 24 * 60 * 60;

starttime_LOS = min(data_LOS(:, 1));
data_LOS(:, 1) = data_LOS(:, 1) - min(data_LOS(:, 1));
data_LOS(:, 1) = data_LOS(:, 1) * 24 * 60 * 60;

%% Filter data by CnR
disp('Filter data by CnR');
cnr_ind = find(cnr_VAD <= -20 | cnr_VAD >= -5);
data_VAD(cnr_ind, 2:3) = NaN;
clear ts_VAD rs_VAD cnr_VAD cnr_ind;

cnr_ind = find(cnr_LOS <= -20 | cnr_LOS >= -5);
data_LOS(cnr_ind, 2:3) = NaN;
clear ts_LOS rs_LOS cnr_LOS cnr_ind;

%% Make time continuous
disp('Make time continuous');
data_LOS_cont = NaN(ceil((max(data_LOS(:, 1))-min(data_LOS(:, 1)))*fs_LOS)+1, length(data_LOS(1, :)));
data_LOS_cont(:, 1) = 0:(1/fs_LOS):(length(data_LOS_cont(:, 1))-1)/fs_LOS;

data_VAD_cont = NaN(ceil((max(data_VAD(:, 1))-min(data_VAD(:, 1)))*fs_VAD)+1, length(data_VAD(1, :)));
data_VAD_cont(:, 1) = 0:(1/fs_VAD):(length(data_VAD_cont(:, 1))-1)/fs_VAD;

for i = 1:length(data_LOS(:, 1))
    data_LOS_cont(round((data_LOS(i, 1)*fs_LOS))+1, :) = data_LOS(i, :);
end

for i = 1:length(data_VAD(:, 1))
    data_VAD_cont(round((data_VAD(i, 1)*fs_VAD))+1, :) = data_VAD(i, :);
end

%% Bin VAD measurements by azimuth angle
binsize = 10;
if (rem(360, binsize~=0))
    disp('Warning fractional bins');
end
rs_binned_VAD = NaN(length(data_VAD_cont(:, 1)), (360/binsize)+2);
rs_binned_VAD(:, 1) = data_VAD_cont(:, 1);

for i = 1:length(data_VAD_cont(:, 1))
    if (~isnan(data_VAD_cont(i, 4)))
        rs_binned_VAD(i, ceil((data_VAD_cont(i, 4)-(binsize/2))/binsize+1)) = data_VAD_cont(i, 2);
    end
end

%% Average over 10 minute intervals
disp('Average over 10 minute intervals');
mean_rs_LOS = NaN(floor(length(data_LOS_cont)/600), 1);
std_rs_LOS = NaN(floor(length(data_LOS_cont)/600), 1);
ts_10min_LOS = NaN(floor(length(data_LOS_cont)/600), 1);
for i = 1:floor(length(data_LOS_cont(:, 1))/600)
    mean_rs_LOS(i) = nanmean(data_LOS_cont((i*600-599):(i*600), 2));
    std_rs_LOS(i) = nanstd(data_LOS_cont((i*600-599):(i*600), 2));
    ts_10min_LOS(i) = data_LOS_cont(i*600-599, 1);
end 

mean_rs_VAD = NaN(floor(length(rs_binned_VAD)/600), length(rs_binned_VAD(1, :))-1);
std_rs_VAD = NaN(floor(length(rs_binned_VAD)/600), length(rs_binned_VAD(1, :))-1);
ts_10min_VAD = NaN(floor(length(rs_binned_VAD)/600), 1);
for i = 1:floor(length(rs_binned_VAD(:, 1))/600)
    mean_rs_VAD(i, :) = nanmean(rs_binned_VAD((i*600-599):(i*600), 2:end));
    std_rs_VAD(i, :) = nanstd(rs_binned_VAD((i*600-599):(i*600), 2:end));
    ts_10min_VAD(i) = rs_binned_VAD(i*600-599, 1);
end 

%% Plot data with outliers
figure();
hold on;
plot(rs_binned_VAD(:, 1), rs_binned_VAD(:, 2), 'gx');

%% Filter outliers
disp('Filter outliers');
for i = 1:600*length(mean_rs_LOS)
    if (data_LOS_cont(i, 2) > (mean_rs_LOS(ceil(i/600)) + 3*std_rs_LOS(ceil(i/600))) || ...
        data_LOS_cont(i, 2) < (mean_rs_LOS(ceil(i/600)) - 3*std_rs_LOS(ceil(i/600))))
        data_LOS_cont(i, 2) = NaN;
    end
end

for i = 1:600*length(mean_rs_VAD(:, 2))
    for j = 2:length(rs_binned_VAD(1, 1:end))
        if (rs_binned_VAD(i, j) > (mean_rs_VAD(ceil(i/600), j-1) + 3*std_rs_VAD(ceil(i/600), j-1)) || ...
            rs_binned_VAD(i, j) < (mean_rs_VAD(ceil(i/600), j-1) - 3*std_rs_VAD(ceil(i/600), j-1)))
            rs_binned_VAD(i, j) = NaN;
            data_VAD_cont(i, 2) = NaN;
        end
    end
end

%% Plot data without outliers
plot(rs_binned_VAD(:, 1), rs_binned_VAD(:, 2), '.b');
plot(ts_10min_VAD, mean_rs_VAD(:, 1)+3*std_rs_VAD(:, 1), '-r');
plot(ts_10min_VAD, mean_rs_VAD(:, 1)-3*std_rs_VAD(:, 1), '-r');
hold off;

%% Fit cosine
disp('Fit cosine');

fitx = (binsize:binsize:(360-binsize))';
type = fittype('a*cosd(x-b)+c');
mycurve = fit(fitx, mean_rs_VAD(120, 2:(end-1))', type, 'StartPoint', [0 0 0]);

%plot(fitx, mean_rs_VAD(120, 2:(end-1)), 'bx');
%hold on;
%plot(fitx, mycurve(fitx), '-r');
%hold off;

%% Plot for comparison
%figure();
%plot(ts_10min_LOS, mean_rs_LOS, '-r');
%hold on;
%plot(data_LOS_cont(:, 1), data_LOS_cont(:, 2), '--b');
%hold off;

%figure();
%plot(ts_10min_VAD, mean_rs_VAD, '-r');
%hold on;
%plot(data_VAD_cont(:, 1), data_VAD_cont(:, 2), '--b');
%hold off;