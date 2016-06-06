%% Task 1
disp('Load Data')
load('WMP_WEnMet_data.mat');
raw_data = readtable('MERRA2_N55.000_E013.125_1992-2016.txt','Delimiter','tab','HeaderLines',26);


fino2_v92 = Fino2.ws92';
fino2_d91 = Fino2.wd91';

last_timestamp = find(Fino2.time==datenum('31-May-2014 23:55:00'))
length(fino2_v92)


% Check for continous time axis
count = 0
for i = 2:length(fino2_v92)
    if Fino2.time(i)-Fino2.time(i-1) >= 1/24/6*1.5
        count = i;
    end
end

fino2_1h_v92 = [];
for i = 1:last_timestamp/6
    fino2_1h_v92(i,1) = datenum('01-Jan-2010 00:00:00')+(i-1)*1/24;
    range_array_v = fino2_v92((i-1)*6+1:i*6,1);
    range_array_dir = fino2_d91((i-1)*6+1:i*6,1);

    fino2_1h_v92(i,2) = nanmean(range_array);
    fino2_1h_v92(i,3) = nanmean(range_array_dir);
end

first_timestamp = find({raw_data{i,1}} == '13.11.2015 08:00');

for 

% timestamps = raw_data.Var1;
% for i = 1:length(timestamps)
%     first_timestamp = strfind('13.11.2015 08:00', timestamps{i});
% end




%fino2_1h_v92 = fino2_1h_v92'
%ino2_1h_v92 = [Fino2.time(1:last_timestamp/6), fino2_1h_v92];
    

% %% Task 3
% wind_sector = [];
% windsOnHeight = [];
% j = 1
% for i = 1:length(fino1_v90)
%     if (fino1_d90(i) >= 240 && fino1_d90(i) <= 285)
%         wind_sector(j,2) = fino1_v90(1,i);
%         wind_sector(j,1) = fino1_d90(1,i);
%         windsOnHeight(j,1) = Fino1.ws33(1,i);
%         windsOnHeight(j,2) = Fino1.ws40(1,i);
%         windsOnHeight(j,3) = Fino1.ws50(1,i);
%         windsOnHeight(j,4) = Fino1.ws60(1,i);
%         windsOnHeight(j,5) = Fino1.ws70(1,i);
%         windsOnHeight(j,6) = Fino1.ws80(1,i);
%         windsOnHeight(j,7) = Fino1.ws90(1,i);
%         windsOnHeight(j,8) = Fino1.ws100(1,i);
%         j= j+1;
%     end
% end

%% ExtraTas