%{

%% Overview
Model: nonlinear_bicycle_model
Test: Step-steer response of vehicle model. The vehicle travels at an
initial velocity and a step impulse of steering angle is applied.

%% Details:
CSAPS based nonlinear tire model.
Fy as function of Fz, SA, IA, P = 8 PSI, SR = 0 (pure cornering)

%% TODO:
- analyze tire data for more accurate relaxation length (has huge effects)
- go through tire model and verify coordinate systems
- validate and improve the vehicle parameters
- IA is modeled as static -2 degrees, this should be pulled from the vehicle parameters
- ensure that inputs into the tire model are bounded (Fz should never go
negative)

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
load("../tire_data/d2704/d2704_7in_csaps.mat");

% open simulink model
open_system('nonlinear_double_track_model');