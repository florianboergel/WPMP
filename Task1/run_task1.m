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
load('data_pp.mat')

% means and 10 min means
disp('Computing 10min means and stddev')
means_interval10 = zeros(round(n/600),20);
for i = 1:n/600
   for j=1:10  
       means_interval10(i,j*2-1) = nanmean(data_pp((i-1)*600+1:i*600,j+1));
       means_interval10(i,j*2) = nanstd(data_pp((i-1)*600+1:i*600,j+1));
   end
end

% Wind directions
u90 = sind(data_pp(:,2));
u33 = sind(data_pp(:,3));
v90 = cosd(data_pp(:,2));
v33 = cosd(data_pp(:,3));
