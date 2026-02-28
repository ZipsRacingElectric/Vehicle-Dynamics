% Codes for MATLAB lessons
% Yaw Gain estamate based on SS Circle + slalom Data from ZR25 Goodyear testing
% Goals: Understnad vehicle response to steering input, extract yaw gain,
% observe how yaw rate changes with speed - under/oversteer

%% Load Data + define time interval

% load in vehicle opject
% addpath vehicle_data ;
% githubFolder = '\vehicle_data\';
% parameterSpreadsheet = strcat(githubFolder,'zr25_data.xlsx');
% ZR25 = vehicle(parameterSpreadsheet);

clear, close all 
clc

% load csv
Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Data Excels\Goodyear data slalom + circles.csv");
Data = Data(Data.timestamps >= 325 & Data.timestamps <= 350, :); % define time range of interest slalom

% load channels
Time = Data.timestamps;
r = -Data.BOSCH_Z_ANGLE_RATE; % deg/s
ay = Data.BOSCH_Y_ACCELERATION; % gs
Speed = Data.SPEED; % Km/hr
Delta = Data.STEERING_ANGLE; % deg


%% Clean and Convert Data

% Clean
r_filt = movmean(r, 10);
Speed_filt = movmean(Speed, 10);
Delta_filt = movmean(Delta, 10);


% Convert Units
r_filt_conv = r_filt*(pi/180); % deg/s -> rad/s
Speed_filt_conv = Speed_filt* (1000/3600);   % km/h -> m/s
Delta_filt_conv = deg2rad(Delta_filt); % deg -> rad

%% Compute Instantaneous Yaw Gain

Yaw_Gain = r_filt_conv./Delta_filt_conv; % compute yaw gain

%% Plots

% sanity data plot
figure()
subplot(3,1,1)
plot(Time, ay)
xlabel('Time')
ylabel('lat accel')
grid on

subplot(3,1,2)
plot(Time, Speed_filt, 'Color', 'r')
xlabel('Time')
ylabel('Speed')
grid on

subplot(3,1,3)
plot(Time, Delta_filt, 'Color','g')
xlabel('Time')
ylabel('Steering ang deg')
grid on



% Yaw Gain Plot 
figure()
scatter(Delta_filt, r_filt)
xlabel('Steering Ang deg')
ylabel('Yaw rate deg/s')
% xlim([0,30])
% ylim([0,30])
grid on

% 
% % Yaw gain vs Speed
% figure ()
% scatter(Speed_filt, Yaw_Gain, 'Color', 'r')
% 
% ylim([-10,10])
% xlabel('Speed')
% ylabel('Yaw Gain')
% grid on 
