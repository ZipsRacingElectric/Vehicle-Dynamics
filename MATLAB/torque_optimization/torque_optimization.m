%{
%% Overview
Test script to validate torque optimization algorithm before implementation
into the control system

%% Details:
This takes an instantaneous snapshot of the vehicle cornering, and a yaw
moment command input. Then it solves for new motor torques that will get us
closer to the target yaw moment

%% TODO:

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

%% Set up cornering situation
% Just using some saved data from the nonlinear_double track model
% response. V = 15 m/s, delta = 20 degrees
load("data/step_steer_sweep_single.mat");

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

% Extract logged data (assuming logsout is used for data logging)
logs = out.logsout;
signal = get(out.logsout, 1); % Replace with actual signal name
vehicle_sensors = logs.get('vehicle_sensors'); % Replace with actual signal name
time = signal.Values.Time;
yaw_ref = signal.Values.Data;
yaw_rate = vehicle_sensors.Values.yaw_rate.Data;
Mz_tv = logs.get('Mz_tv').Values.Data;

% Plot transient response of run
figure;
plot(time, yaw_ref, 'r--', 'LineWidth', 1.5); % Yaw reference (dashed red line)
hold on;
plot(time, yaw_rate, 'b-', 'LineWidth', 1.5); % Yaw rate (solid blue line)

% Labels and title
xlabel('Time (s)');
ylabel('Yaw (deg/s)');
title('Yaw Reference vs. Yaw Rate');
legend('Yaw Reference', 'Yaw Rate');
grid on;
hold off;

%% Grab instanteous data for testing torque optimization
index = 55;
t = time(index);
apps = 0;

%% Instantaneous Yaw Moment and Yaw Moment Target
% 1st objective we will try to achieve
% We will take the derivative of the measured yaw rate sensor
time_step = time(index) - time(index - 1);
[yaw_moment, ~] = get_yaw_moment(yaw_rate(index), yaw_rate(index - 1), time_step, zr25.yaw_polar_inertia);
yaw_moment_tv = Mz_tv(index);

%% Instantaneous Long. Acceleration and Long Acceleration Target
% 2nd object we will try to achieve
x_dot = vehicle_sensors.Values.body_acceleration_x.Data;
x_dot_ref = get_x_accel_target(apps);

%% Calculate Limitations of the Powertrain
% The optimization needs to be constrained around these limits




% Calculate yaw moment from derivative of yaw rate sensor
function [yaw_moment, yaw_acceleration]  = get_yaw_moment(yaw_rate, yaw_rate_prior, time_step, inertia)
    yaw_acceleration = yaw_rate - yaw_rate_prior / time_step; % rad/s^2
    yaw_moment = yaw_acceleration * inertia;
end


% Calculate longitudinal acceleration target from driver sensors
function x_dot_ref  = get_x_accel_target(apps)
    % placeholder for some more complicated algorithm in the future
    x_dot_max = 2 * 9.81; % m/s^2
    x_dot_ref =  x_dot_max * apps / 100;
end