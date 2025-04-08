%{
%% Overview
This script runs a simulink model of a simple PID system used to verify the
hardware-in-the-loop (HIL) system with the VCU is functional and works as
expected, before testing torque vectoring control algorithms. It interfaces
with the hardware during simulation by sending and recieving CAN messages
with the device under test.

%% Details:
- Not real-time HIL
- im running python in a virtual enviornment in the repository directory
/venv/ - set up your virtual python enviornment as necessay for your system
- you may need to make the script executable with additional system
commands for your given system

%% TODO:
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

%% Select Simulation
MODEL = 'sine_sweep_model';
%MODEL = 'step_steer_model';

%% Load Model Data
% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path
% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");
% create simulink parameters
zr25.create_simulink_parameters();
% load tire model
% Note: this data is fixed at P = 8 psi, IA = 2 deg
load("../tire_data/d2704/d2704_7in_rbf_full_model.mat");
% load aero data
load("../vehicle_data/zr25_aero_data.mat");
% load steering data
load("../vehicle_data/zr25_steering_data.mat");

test = false;

%% Test that the python can interface is working
if (test)
    input_data = -0.2354;
    plant_data = 10.6;
    time_step  = 300;
    vector = [input_data, plant_data, time_step];
    output = send_can_message([1000 -8.560702e+02 10]);
    output = send_can_message([1000 -8.560702e+02 10]);
    output = send_can_message([1000 -8.560702e+02 10]);
    output = send_can_message([1000 -8.560702e+02 10]);
    output = send_can_message([1000 -8.560702e+02 10]);
    output = send_can_message('exit');
else
    open_system(MODEL);
end