%{
%% Overview
Model: nonlinear_double_track_model
Test: Step-steer response of vehicle model. The vehicle travels at an
initial velocity and a step impulse of steering angle is applied.

This model will serve as the basis for the LPV plant model used for torque
vectoring control.

%% Details:
RBF-based nonlinear tire model.
Fy as function of SA, FZ, IA, P = 8 PSI, SR = 0 (pure cornering)
Mz as a function of SA, FZ, IA, P = 8 PSI, SR = 0 (pure cornering)
- Positive steering wheel input is turning the wheel to the right.

%% TODO:
- electric motor model (torque request input, inertia, wheel speed, and
slip ratio output)
- validate the vehicle parameters
- velocity input for Mz should be the tire longitudinal velocity and not the
longitudinal velocity for the vehicle
- adjust rear Fy steer compliance based on toe link location & trail (under or over steer
induced compliance, depends on design)
- linear scaling factor for tire-road mu correction should be changed to a
non-linear gain based on a constant radius test of the car (see tire
correction document in the tire_data/tire model research/ folder)
- add Slip Ratio calculation
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

% % load vehicle data
% zr25 = vehicle("../vehicle_data/zr25_data.xlsx");
% 
% % create simulink parameters
% zr25.create_simulink_parameters();
% 
% % load tire model
% % Note: this data is fixed at P = 8 psi, IA = 2 deg
% load("../tire_data/d2704/d2704_7in_rbf_full_model.mat");
% 
% % load aero data
% load("../vehicle_data/zr25_aero_data.mat");
% 
% % load steering data
% load("../vehicle_data/zr25_steering_data.mat");

% open simulink model
open_system('pid_testing_model');

% Set simulation parameters
sim('pid_testing_model');  % Runs the simulation

load('inputData.mat');
load('outputData.mat');

inData=inputData.'; %transpose
outData=outputData.';

inData(:,1)=[];
outData(:,1)=[];

% Now save the workspace variables
save('dataFile.mat', 'inData', 'outData');

load dataFile
Z = iddata(outData,inData,0.05);
Z.InputName = 'step input';
Z.InputUnit = 'rad';

Z.OutputName = {'ouput'};
Z.OutputUnit = {'rad'};

t = Z.SamplingInstants;

subplot(2,1,1);
plot(t,Z.y), ylabel('Output (rad)')
title('Logged Input-Output Data')
% axis([0 50 -1.2 1.2])

subplot(2,1,2);
plot(t,Z.u), ylabel('Input (rad)')
% axis([0 50 -1.2 1.2])
xlabel('Time (seconds)')

opt = ssestOptions('Focus','simulation');
syslin2 = ssest(Z, 2, 'DisturbanceModel', 'none', opt);
syslin3 = ssest(Z, 3, 'DisturbanceModel', 'none', opt);
syslin4 = ssest(Z, 4, 'DisturbanceModel', 'none', opt);
syslin5 = ssest(Z, 5, 'DisturbanceModel', 'none', opt);
% choose the best option

syslin.InputName = Z.InputName;
syslin.OutputName = Z.OutputName; % reconcile names to facilitate comparison
clf
compare(Z, syslin2, syslin3, syslin4, syslin5)

% Extract state-space matrices from the idss object
A = syslin5.A;  % Extract state matrix
B = syslin5.B;  % Extract input matrix
C = syslin5.C;  % Extract output matrix
D = syslin5.D;  % Extract direct transmission matrix

% create the transfer function from the state-space representation
[num, den] = ss2tf(A, B, C, D)