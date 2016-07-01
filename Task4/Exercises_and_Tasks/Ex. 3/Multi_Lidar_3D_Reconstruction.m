close all; clc;
      
% Measurement geometry [x, y, z] of the measurement point and the three
% Lidars

staring_point    =   [ 46.8,     0, 90]';
lidar_positions  =   [-0.94,-34.59, 0.60;
                      -1.01, 49.58, 2.09;
                      78.69,  5.40, -0.74]';
         
lidar_range=sqrt((repmat(staring_point(1),1,3)-lidar_positions(1,:)).^2+...
           (repmat(staring_point(2),1,3)-lidar_positions(2,:)).^2);
       
% Here you have to calculate the azimuth and elevation scanning angles of each of
% the lidars:
for i=1:3
    distOnGround = (staring_point(1) - lidar_positions(i,1))
   ele(i)= acos();
end


% Reading the line-of-sight velocity data

vlos1 = csvread('R2D1.csv',1,0);
vlos2 = csvread('R2D2.csv',1,0);
vlos3 = csvread('R2D3.csv',1,0);

vlos1(:,1)=vlos1(:,1)/1000;     % First column is time in milliseconds
vlos2(:,1)=vlos2(:,1)/1000;     % First column is time in milliseconds
vlos3(:,1)=vlos3(:,1)/1000;     % First column is time in milliseconds

fs = 97.65625;  % Sampling frequency
dt = 1/fs;      % Time step
starttime = dt;
endtime   = min([vlos1(end,1) vlos2(end,1) vlos3(end,1)]);
timeline = (starttime:dt:endtime)';

% The line-of-sights speeds are interpolated to the same timeline to be
% perfectly synchronised. Continue working with these vectors:

vlos1_interp = interp1(vlos1(:,1),vlos1(:,2),timeline);
vlos2_interp = interp1(vlos2(:,1),vlos2(:,2),timeline);
vlos3_interp = interp1(vlos3(:,1),vlos3(:,2),timeline);

% Conversion of the three line-of-sight speeds to u and v

% matrix = ...
% u = ...
% v = ...
% w = ...

% In this part the calculated u, v and w components are aligned with the
% 'main' wind direction, such that both v and w have a mean of 0 m/s:

theta = atan(mean(v)/mean(u));
delta = -atan(mean(w)/sqrt(mean(u)^2+mean(v)^2));
alignmatrix1=Rotate_Matrix(delta*180/pi,'y');
alignmatrix2=Rotate_Matrix(theta*180/pi,'z');
V_aligned=alignmatrix2*alignmatrix1*V'; V_aligned=V_aligned';
statistics = [ mean(V_aligned);
                std(V_aligned)];
            
% Plotting of results

colors='bgr';
figure
hold on
 for i=1:3
    scatter3(lidar_positions(1,i),...
          lidar_positions(2,i),...
          lidar_positions(3,i),'k','linewidth',3)
      textstring=sprintf('%s%d','R2D',i);
      text(lidar_positions(1,i),...
          lidar_positions(2,i)-5,textstring,'fontsize',12);
          formatstring=sprintf('%s%s',colors(i),'--');
    plot3([lidar_positions(1,i),staring_point(1)*1.3-0.3*lidar_positions(1,i)],...
          [lidar_positions(2,i),staring_point(2)*1.3-0.3*lidar_positions(2,i)],...
          [lidar_positions(3,i),staring_point(3)*1.3-0.3*lidar_positions(3,i)],...
          formatstring,'linewidth',2);
 end
 grid on
 axis equal
 xlabel('x [m]','fontsize',14)
 ylabel('y [m]','fontsize',14)
 zlabel('z [m]','fontsize',14)
 set(gca,'fontsize',12)
 xlim([-5 80]); ylim([-40 50]); zlim([-1 120]);
 view([-35 20])