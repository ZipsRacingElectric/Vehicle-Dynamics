%{
%% Overview
This script runs a simulink model of a simple PID system used to verify the
hardware-in-the-loop (HIL) system with the VCU is functional and works as
expected, before testing torque vectoring control algorithms. It interfaced
with the hardware during simulation by sending and recieving CAN messages
with the device under test.

%% Details:
- Not real-time HIL
- am running python in a virtual enviornment in the repository directory
/venv/ - set up your virtual python enviornment as necessay for your system
- you may need to make the script executable with additional system
commands for your given system

%% TODO:
%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

test = false;

%% Test that the python can interface is working
if (test)
    input_data = -0.2354;
    plant_data = 10.6;
    time_step  = 300;
    vector = [input_data, plant_data, time_step];
    
    % This function is on the path in send_can_message.m
    output = send_can_message(vector);

else
    open_system('hil_simulink');
end