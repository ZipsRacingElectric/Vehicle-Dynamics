%{
%% Overview
This script runs a simulink model of a simple PID system used to verify the
hardware-in-the-loop (HIL) system with the VCU is functional and works as
expected, before testing torque vectoring control algorithms.

%% Details:
- MATLAB must be setup to work with a python enviornment. See more at
https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html
- Note that MATLAB does not have support for the latest python 3 version. I
am running the latest python version supported by R2023b (pytohn 3.11) in a
virtual enviornment in the repository directory /venv/ - set up your virtual
python enviornment as necessay for your system

%% TODO:
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% Setup the python interpreter to use the python venv:
pyenv('Version', './venv/bin/python')

% Add the current MATLAB working directory to Python's sys.path
currentDir = pwd;
if count(py.sys.path, currentDir) == 0
    insert(py.sys.path, int32(0), currentDir);
end

try
    py.importlib.import_module('can');
    disp('python-can imported successfully');
catch ME
    disp(ME.message);
end

% Attempt to import the module
try
    py.importlib.import_module('can_module')
catch e
    disp('Error importing the module:');
    disp(e.message);
end

py.print('Hello from Python!')

input_data = int32(10);
plant_data = int32(20);
time_step  = int32(30);

actuating_signal = pyrunfile('can_module.py', 'actuating_signal', ...
    'input_data', input_data, 'plant_data', plant_data, 'time_step', time_step);
disp(actuating_signal);
%open_system('hil_simulink');