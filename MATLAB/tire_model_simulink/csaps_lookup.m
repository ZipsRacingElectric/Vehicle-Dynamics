%{
Model: csaps tire splines implemented as a lookup table in simulink
Test: basic sweeps to verify spline data

Details:
Save a .mat data file containing the lookup data at the filepath
indicated below. Csaps spline lookup data should be fitted using the .m script located
in the tire_data/d2704/ directory, or using information availible on the
FSAE TTC Forum.

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load tire data
load("../tire_data/d2704/d2704_7in_csaps.mat");

% open simulink model
open_system('csaps_lookup_model');
