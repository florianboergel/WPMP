%% Spectra_Data
%  Script that converts raw spectra from the SpinnerLidar to the
%  line-of-sight wind speeds

close all; clc;

%% Reading the line-of-sight data of the SpinnerLidar

spinnerlidar_data    = dlmread('SpinnerLidar_Data_1s.txt');
spinnerlidar_spectra = dlmread('SpinnerLidar_Spectra_1s.txt');

index       = spinnerlidar_data(:,1);  % Lidar measurement index
vlos      	= spinnerlidar_data(:,3);  % Line-of-sight measurement
sx          = spinnerlidar_data(:,7);  % Laser pointing unit vector x-component
sy         	= spinnerlidar_data(:,8);  % Laser pointing unit vector y-component
focus       = spinnerlidar_data(:,9);  % Focus distance of Lidar

sz          = sqrt(1-sx.^2-sy.^2);
x           =  sz.*focus;
y           = -sy.*focus;
z           =  sx.*focus;
               
%% Reading the spectrum data of the SpinnerLidar

bins     = 256;             % amount of bins in each spectrum
bandwith = 25e6;            % frequency bandwith of the spectra
lambda  = 1560e-9;          % wavelength of the laser light

%% Step 1: noise cancelling by discarding all bins below mean*1.1 per spectra
for pos = 1:312
    
    spinner_noiseCancelled(pos,:) = spinnerlidar_spectra(pos,:)-1.5*mean(spinnerlidar_spectra(pos,:));
    spinner_noiseCancelled(spinner_noiseCancelled < 0) = 0;
%    figures created once in directory
%    figure('visible','off')
%    plot(spinner_noiseCancelled(pos,:));
%    saveas(gcf,strcat('figures/spectra_noisecancels_normed_',num2str(pos) ,'.jpg'));
end;

%% Step 2
% Calculate the frequency resolution of the spectra, the Doppler frequency of each
% bin and from that the wind speed corresponding to each bin:
c = 3*10^8;
f_0 = c/lambda;
lambda_0 = c/f_0;
for bin=1:256
    f_d(bin,1) = (bin-1)/bins*bandwith;
    v(bin,1) = f_d(bin,1) *lambda_0; %divided by 2?
end;

%% Finding the spectra peaks and visualize
%%Step 3: apply centroid method
for pos = 1:312 
    spinnerSum = sum(spinner_noiseCancelled(pos,:))
    spinner_noiseCancelled_normed(pos,:) = spinner_noiseCancelled(pos,:)/spinnerSum;
end;

for pos = 1:312
    f_peak(pos,1) = sum(spinner_noiseCancelled_normed(pos,:)*f_d(:))/sum(spinner_noiseCancelled_normed(pos,:));
end;

%detect centroid failures
for pos = 1:312
    maxIndex= find(spinner_noiseCancelled_normed(pos,:) == max(spinner_noiseCancelled_normed(pos,:))); 
    centroidIndex= f_peak(pos)/bandwith*bins;
    failures(pos,1) = abs(centroidIndex-maxIndex);
end;
figure();
hold on;
x_range = 1:length(failures)
% failures_mean = nanmean(failures);
% failures_corrected = failures - failures_mean;
% failures(failures_corrected < 0) = 0;
plot(x_range,failures);
for k=1:length(failures)
      if failures(k) > 5
      text(x_range(k),failures(k),num2str(x_range(k)))
      end
end
%[pks,locs] = findpeaks(failures);
%text(locs+.02,pks,num2str((1:numel(pks))'));
xlabel('Spectra number','fontsize',10)
ylabel('Differences of bin indices for max and centroid functions','fontsize',10)
saveas(gcf,strcat('figures/failures.png'));
hold off;

%%Step 4: Calc centroid velocity for correlation with lidar measured speeds
vlos_centroid(:,1) = f_peak(:,1) *lambda_0; %divided by 2?

%% Plotting
figure();
plot(spinnerlidar_spectra(1,:))
xlabel('Frequency')
ylabel('Collection efficiency')
title('Spectrum without filter')
saveas(gcf,strcat('figures/spectra_nofilter.png'));

figure();
plot(spinner_noiseCancelled(1,:))
title('Spectrum with filter')
xlabel('Frequency')
ylabel('Collection efficiency')
saveas(gcf,strcat('figures/spectra_noisecanceld.png'));

check = vlos-vlos_centroid;
Correlation(vlos,vlos_centroid);

Rosette_Scan_Plot(y,z,vlos_centroid,...
        'coloraxis',[5 10],...
        'rosettebackground',y,z,...
        'meshgridvector',(-1:0.1:1)*unique(focus)/2);

Rosette_Scan_Plot(y,z,vlos,...
        'coloraxis',[5 10],...
        'rosettebackground',y,z,...
        'meshgridvector',(-1:0.1:1)*unique(focus)/2);
    
  