%{
Model: Basic Load transfer with aero forces
Test: FSAE acceleration event

Details:

More accurate acceleration event than OptimumLap can do. Adds load
transfer, aero forces, and rollout distance per the event rules (rollout
makes a huge difference).

TODO:
- add nonlinear tire model for tire load sensitivity
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load vehicle simulink parameters
run("../vehicle_data/zr25.m");

% load tire model

% open simulink model
open_system('acceleration_event_simulink.slx');