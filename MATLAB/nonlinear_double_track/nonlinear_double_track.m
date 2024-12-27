%{

%% Overview
Model: nonlinear_bicycle_model
Test: Step-steer response of vehicle model. The vehicle travels at an
initial velocity and a step impulse of steering angle is applied.

This model will serve as the basis for the LPV plant model used for torque
vectoring control.

%% Details:
CSAPS based nonlinear tire model.
Fy as function of SA, FZ, IA, P = 8 PSI, SR = 0 (pure cornering)
Mz as a function of SA, FZ, IA, P = 8 PSI, SR = 0 (pure cornering)
- Positive steering wheel input is turning the wheel to the right.

%% TODO:
- go through tire model and verify coordinate systems
- validate and improve the vehicle parameters
- velocity input for Mz should be the tire longitudinal velocity and not the
longitudinal velocity for the vehicle
- basic front compliance (Mz induced steer compliance)
- basic rear compliance (Fy induced steer compliance (toe compliance))
- investigate lateral load transfer distribution and see if it can be
easilly improved using calculated values from spring, rollbar, & damper
design (currently just 50:50)
- linear scaling factor for tire-road mu correction should be changed to a
non-linear gain based on a constant radius test of the car (see tire
correction document in the tire_data/tire model research/ folder)

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

% load tire model
tire_fy = load("../tire_data/d2704/d2704_7in_csaps.mat");
tire_mz = load("../tire_data/d2704/d2704_7in_csaps_mz.mat");

% open simulink model
open_system('nonlinear_double_track_model');