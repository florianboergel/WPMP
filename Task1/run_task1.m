%% preprocess data
disp('Preprocessing Data')
disp('Loading Data ...')
raw_data = readtable('1301.txt','Delimiter','tab');
disp('Writing timestamp ...')
time_stamp = raw_data{:, {'Time'}};
disp('Writing raw_data')
raw_data = raw_data{:, {'d90', 'd33', 'u100', 'u90', 'u80', ...
'u70', 'u60', 'u50', 'u40', 'u33'}};
disp('Converting time with datenum ...')
t = (datenum(time_stamp, 'yyyy-mm-dd HH:MM:SS'));
clear time_stamp;
t = round(t*24*3600);
disp('Marking invalid data ...')
%raw_data(raw_data==-999) = NaN(size(raw_data(raw_data==-999)));
raw_data(raw_data==-999) = NaN;
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

%% Evaluation
load('data_pp.mat');

% means and 10 min means
disp('Computing 10min means and stddev');
means_interval10 = zeros(round(n/600),20);
for i = 1:n/600
   %treat wind vanes
   radians=data_pp((i-1)*600+1:i*600,2)/180*pi;
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

start30thJan = 29*24*6;
plot(1:24*6, means_interval10(start30thJan+1:start30thJan+24*6,5),1:24*6,means_interval10(start30thJan+1:start30thJan+24*6,5)-means_interval10(start30thJan+1:start30thJan+24*6,6)/2,1:24*6,means_interval10(start30thJan+1:start30thJan+24*6,5) + means_interval10(start30thJan+1:start30thJan+24*6,6)/2);
save('meansAndStddev.mat', 'means_interval10');