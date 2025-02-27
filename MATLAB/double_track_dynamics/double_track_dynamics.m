%{
%% Overview
Model: double_track_dynamics
Test: Sweep some tire forces to validate model

This model can be used to replace the more complicated double track block
provided by MATLAB. It also uses coordinate systems in alignment with the
rest of the torque vectoring design

%% Details:
This block solves 3 equations of motion for each degree of freedom.

- 3 degrees of freedom (X, Y, yaw axis)
- fixed inertial frame outputs (position_x, position_y, yaw_angle,
yaw_rate)
- body reference frame outputs (body_velocity_x, body_velocity_y)

%% TODO:
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

% create simulink parameters
zr25.create_simulink_parameters();

% open simulink model
open_system('double_track_dynamics_model');