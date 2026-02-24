%% Understeer Gradient Calculator
%The understeer gradient is a key metric in vehicle dynamics that describes how a car’s steering response changes with lateral acceleration. It quantifies whether a vehicle tends to understeer (requires more steering angle as cornering force increases), oversteer (requires less steering angle), or is neutral steer.
%A positive understeer gradient means the car needs progressively more steering input to hold higher lateral acceleration → typical for production cars (stable, predictable).
%A negative gradient indicates oversteer → the car needs less steering input as cornering builds, which can become unstable at the limit.
%Zero gradient is neutral steer → steering angle vs. lateral acceleration is linear.
%This metric is useful because it:
% Provides a simple, quantitative way to describe handling balance.
%Helps engineers compare setup changes (springs, bars, tire pressures, aero balance).
%Links test data (steering, lateral accel, speed) to fundamental handling theory.
%Serves as a tuning target — FSAE teams often aim for a mild understeer gradient for stability while still keeping sharp turn-in.

%Constant speed radius
%(SteerAngle/LatAccel)-(Wheelbase/velocity^2)

clc, clear, close all
% Filepath to CSV
filepath = "C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Data Excels\Understeer Gradient GY 02_19_26 ZR25.csv";

% Load table
Data = readtable(filepath);
Data = Data(Data.timestamps >= 130 & Data.timestamps <= 155, :); % define time range of interest
time = Data.timestamps; % s

% load in vehicle object 🏎️
addpath vehicle_data ;
githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr25_data.xlsx');
ZR25 = vehicle(parameterSpreadsheet);

L = ZR25.wheelbase;
steeringratio = ZR25.steering_ratio;
g = 9.81;

%% Extract steering angle and lateral acceleration
delta_raw = Data.STEERING_ANGLE;  % degrees 
U_raw  = Data.SPEED;                    % km/h
r_raw  = Data.BOSCH_Z_ANGLE_RATE;       % deg/s
ay_raw = Data.BOSCH_Y_ACCELERATION;     % g

% Unit Conversions 
U_raw_con  = U_raw * (1000/3600);   % km/h -> m/s
r_raw_con  = r_raw * (pi/180);      % deg/s -> rad/s
ay_raw_con = ay_raw * 9.81;   % g -> m/s^2
delta_raw_con = deg2rad((delta_raw)/steeringratio); %deg - rad

%  Low-Pass Butterworth Filter
Fs = 100;
fc = 8;                         % Cutoff frequency (Hz)
[bfilt,afilt] = butter(2, fc/(Fs/2));    % 2nd-order Butterworth

% Apply Zero-Phase Filtering 
U  = filtfilt(bfilt,afilt,U_raw_con);
r  = filtfilt(bfilt,afilt,r_raw_con);
ay = filtfilt(bfilt,afilt,ay_raw_con);
delta = filtfilt(bfilt, afilt, delta_raw_con);


%%  STEADY-STATE MASK
r_dot = [0; diff(r)] * Fs;
delta_rate = [0; diff(delta)] * Fs;
ay_dot = [0; diff(ay)] * Fs;

mask = ...
    abs(r_dot) < 0.05 & ...
    abs(delta_rate) < deg2rad(20) & ...
    abs(ay_dot) < 0.4 & ...
    U > 4 & ...
    abs(r) > 0.15 & ...      
    abs(ay/g) > 0.3;          

% apply mask
U_ss     = U(mask);
r_ss     = r(mask);
delta_ss = delta(mask);
ay_ss    = ay(mask);
time_ss  = time(mask);


%%  POINTWISE Ku 
K_point = (g ./ U_ss) .* (delta_ss ./ r_ss - L ./ U_ss);  % rad/g
K_deg   = rad2deg(K_point);

%% GLOBAL ku
A = [(L ./ U_ss) .* r_ss, (U_ss / g) .* r_ss];
theta = A \ delta_ss;

Ku_global = rad2deg(theta(2));

fprintf('\nGlobal Understeer Gradient = %.3f deg/g\n', Ku_global);


%% PLOT #1 (PRIMARY — driver meaningful) 
figure;
scatter(ay_ss/g, rad2deg(delta_ss), 25, 'filled')
xlabel('Lateral Acceleration [g]')
ylabel('Road Wheel Steering Angle [deg]')
title('Steering vs Lateral Acceleration')
grid on


figure;
scatter(ay_ss/g, rad2deg(delta_ss), 25, 'filled')
xlabel('Lateral Acceleration [g]')
ylabel('wheel steering ang deg')
title('Pointwise Understeer Gradient')
grid on


figure;
scatter(time_ss, ay_ss/g, 20, 'filled')
xlabel('Time [s]')
ylabel('Lateral Acceleration [g]')
title('Filtered Steady-State Segment')
grid on
