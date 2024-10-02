%{
Model: csaps tire splines implemented in simulink
Test: basic sweeps to verify spline data

Details:
Save a .mat data file containing the csaps splines at the filepath
indicated below. Csaps splines should be fitted using the .m script located
in the tire_data/d2704/ directory, or using information availible on the
FSAE TTC Forum.

Matlab csaps spline functions fit splines to 3D surface data, with 1 output
and 2 inputs. This means that for each surface, usually tire pressure and
IA are fixed, and slip angle and Fz are the input parameters. Csaps fits
pretty well compared to pacejeka, is quicker to evaluate, and can be
inversed easilly to use Fy or Fx as an input parameter.

This simulink model uses a single csaps surface to provide tire data,
however with a kinematic model it may not be that useful since you will
need 3 inputs (Fz, slip angle, IA). Building off this, a model which
interpolates output values between multiple csaps surfaces built off of fixed
values of IA is the next step.
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% load tire data
run("../tire_data/d2704/d2704_csaps.mat");

% open simulink model
% add an interpreted matlab function block that calls the following:
% fnval(MY_SURFACE_NAME,{u(1),u(2),u(3)})
open_system('csaps_tire_model');
