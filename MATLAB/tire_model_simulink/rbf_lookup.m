%{
Model: lookup tire model generated via rbf interpolation, implemented in
simulink. This model includes both Fy and Fx effects
Test: basic sweeps to verify spline data

Details:
- Save a .mat data file containing the csaps splines at the filepath
indicated below. Csaps splines should be fitted using the .m script located
in the tire_data/d2704/ directory
- This model uses a simplified lookup table to validate the RBF method. The
lookup table assumes 8 psi and 2 deg IA, with SA, SL, Fz inputs. The
outputs are Fy, Fx, Mz
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load tire data
load("../tire_data/d2704/d2704_7in_rbf_full_model.mat");

% open simulink model
% add an interpreted matlab function block that calls the following:
% fnval(MY_SURFACE_NAME,{u(1),u(2),u(3)})
open_system('rbf_lookup_model');
