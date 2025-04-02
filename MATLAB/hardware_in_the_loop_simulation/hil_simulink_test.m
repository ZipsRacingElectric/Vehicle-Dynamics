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

test = true;

%% Test that the python can interface is working
if (test)
    input_data = int32(10);
    plant_data = int32(20);
    time_step  = int32(30);
    
    % Construct the command string.
    % Make sure to use the correct path to your Python interpreter in the virtual env.
    command = sprintf('./venv/bin/python can_module.py %d %d %d', input_data, plant_data, time_step);
    
    % Call the Python script using system()
    [status, cmdout] = system(command);
    
    % Check if the command executed successfully.
    if status == 0
        % Assume the final line of output contains the actuating signal.
        % Split the output by newline and take the last non-empty line.
        lines = strsplit(strtrim(cmdout), '\n');
        actuating_signal_str = lines{end};
        actuating_signal = str2double(actuating_signal_str);
        disp(cmdout);
        fprintf('Actuating Signal: %f\n', actuating_signal);
    else
        disp('Error running Python script:');
        disp(cmdout);
    end

else
    open_system('hil_simulink');
end