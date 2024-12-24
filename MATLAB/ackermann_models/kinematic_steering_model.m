%{

%% Overview
Model: kinematic_steering

%% Details:
This implements the matlab block that models the steering response as a
function of the rack and pinion and steering geometry. It appears not to
take in caster angle of the kingpin axis, so yeah.

Evaluation as a sub-block in a larger dynamic model.

- The input is at the steering rack input, where the steering sensor
measures. U-joint non-linearity isnt measured by the sensor, so it is not
included.
- Positive steering input is turning the wheel to the right.
- A positive angle at the tire means it is turning right.
- Negative static toe is toe out, and it is per tire.

%% TODO:
- output of steering wheel motion as function of u-joint design
- values are from ZR25_SuspensionForces.xlsx, need to confirm their validity

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

%% Calculate some relevant suspension parameters
%{
ubj = [0.804805, 0.554950, 0.340000]; % upper ball joint X, Y, Z in [m]
lbj = [0.817649, 0.577022, 0.130000]; % lower ball joint
tie = [0.891342, 0.578549, 0.170000];
chassis_tie = [0.870000, 0.192088, 0.170000];

kingpin_vector = ubj - lbj % vector of kingpin angle
u_kingpin = kingpin_vector ./ norm(kingpin_vector); % unit vector of kingpin
tie_ubj = tie - ubj; % vector from any point on the kingpin line to the tie rod point

steering_arm_length = Simulink.Parameter(norm(cross(tie_ubj, u_kingpin))); % perpendicular length between tie rod point and kingpin angle
rack_length = Simulink.Parameter(abs(2 * chassis_tie(2))); % length of steering rack
front_tie_rod_length = Simulink.Parameter(norm(tie - chassis_tie)); % length of front tie rod

kingpin_center = (ubj + lbj) ./ 2;
rack_to_axis_distance = Simulink.Parameter(abs(kingpin_center(1) - chassis_tie(1))); % kind of an estimate since the definition assumes a 2d steering system without 3d kingpin axis
steering_pinon_radius = Simulink.Parameter(0.0229); % [m], from ZR20, I assume we havent changed it since
%}

% open simulink model
open_system('kinematic_steering');