%{
%% Overview
Model: aero_model
Test: Sine wave input for velocity into aero lookup tables

%% Details:
This model just verifies the lookup tables for aero data work correctly.


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

% load aero data
load("../vehicle_data/zr25_aero_data.mat");

% open simulink model
open_system('aero_model');