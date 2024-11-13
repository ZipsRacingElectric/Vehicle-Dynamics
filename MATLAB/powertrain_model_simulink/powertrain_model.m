%{
Model: powertrain_model_simulink
Test: N/A

Details:

Full powertrain model including battery pack, inverter, electric motors,
and drivetrain. This can be integrated into a larger vehicle simulation.

TODO:
- 

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load vehicle simulink parameters
run("../vehicle_data/zr25.m");

% load battery pack simulink parameters
run("../battery_data/zr25_battery_pack.m");

% open simulink model
open_system('powertrain_model_simulink.slx');