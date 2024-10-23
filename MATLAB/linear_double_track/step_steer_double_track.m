% Model: linear_double_track_model
% Test: Step-steer response of vehicle model
%{
Details:
From 2024 Michigan Skidpad: 1.28 rad/s @ 20.52 m/s is the time to beat
We are sweeping step steer angles in radians to find the passive vehicle
response.
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load vehicle simulink parameters
run("../vehicle_data/zr25.m");

% open simulink model
open_system('linear_double_track_model');

% steering wheel angle min max sweep values in radians
% positive is steering right?
steering_angle_min = deg2rad(1);
steering_angle_max = deg2rad(180);
num_sweeps = 10;

% parameter sweep of steering angle
for k = 1:num_sweeps
    angle = steering_angle_min + k*((steering_angle_max - steering_angle_min) / (num_sweeps - 1));
    set_param([bdroot,'/steering_kinematics/steering_wheel_angle'], 'Value', num2str(angle));
    simOut(k) = sim(gcs);
    disp(['Completed ', num2str(k),' of ',num2str(num_sweeps),' simulations.']);
end


% Plot step-steer results
for k = 1:num_sweeps
    logsout = simOut(k).logsout;
    Output = logsout.getElement('YawRate');
    y_out = Output.Values.Data;
    x_out = Output.Values.Time;
    plot(x_out, y_out);
    hold on;
end
grid on;
xlabel('Seconds');
ylabel('Yaw-Rate (rad/s)');
legend('Steering wheel 1-180 deg, 20 increments',location = 'northwest');
title('Yaw-rate Step Response');
