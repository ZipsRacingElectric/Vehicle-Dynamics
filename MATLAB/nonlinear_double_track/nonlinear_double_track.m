%{
Model: nonlinear_bicycle_model
Test: Step-steer response of vehicle model

Details:

CSAPS based nonlinear tire model.
Fy as function of Fz, SA, IA, P = 8 PSI, SR = 0 (pure cornering)

TODO:
- analyze tire data for more accurate relaxation length (has huge effects)
- go through tire model and verify coordinate systems
- validate and improve the vehicle parameters

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load vehicle simulink parameters
run("../vehicle_data/zr25.m");

% load tire model
load("../tire_data/d2704/d2704_7in_csaps.mat");

% open simulink model
open_system('nonlinear_bicycle_model');