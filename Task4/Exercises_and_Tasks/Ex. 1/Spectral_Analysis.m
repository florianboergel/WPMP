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

% Calculate the frequency resolution of the spectra, the Doppler frequency of each
% bin and from that the wind speed corresponding to each bin:

% f = ...
% v = ...


%% Finding the spectra peaks




%% Plotting

Correlation(vlos,vlos);

Rosette_Scan_Plot(y,z,vlos,...
        'coloraxis',[5 10],...
        'rosettebackground',y,z,...
        'meshgridvector',(-1:0.1:1)*unique(focus)/2);