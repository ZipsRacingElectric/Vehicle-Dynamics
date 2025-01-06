%{

%% Overview
Model: tire_normal_loads_model
Test: sine sweep of vehicle forces

This model calculates tire normal loads from load transfer, aero forces, 
and self-aligning moment, based on applied vehicle forces, vehicle geometry,
and total lateral load transfer distribution (TLLTD). This model does not 
consider the delay in tire force transfer due to the supsension spring 
and damper system (it is not dynamic, just final steady state values).

%% Details:
- (+) Fx is forward acceleration
- (+) Fy is acceleration when turning right
- Tire forces saturate so Fz never becomes negative

%% TODO:
- output to indicate when the vehicle is lifting a tire
- tire self aligning moment Mx
- aero forces

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
%{
tire_fy = load("../tire_data/d2704/d2704_7in_csaps.mat");
tire_mz = load("../tire_data/d2704/d2704_7in_csaps_mz.mat");
%}

% open simulink model
open_system('tire_normal_loads_model');