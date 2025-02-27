%{

%% Overview
Model: nonlinear_bicycle_model
Test: Sweep step-steer test at various operating conditions of the vehicle.
The starting velocity and the steering angle are swept. 

%% Details:
This uses the Simulation Manager to batch run multiple simulations to get 
an idea of the full response of the vehicle.

%% TODO:
- identify runs which the car is not stable
- Plot the input conditions as inputs and view the region of stability of
the car
- in transient plot, label the legends as combinations of input states
(instead of just run number which is not useful)

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

%% Load Model
% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");
% create simulink parameters
zr25.create_simulink_parameters();
% load tire model
% Note: this data is fixed at P = 8 psi, IA = 2 deg
load("../tire_data/d2704/d2704_7in_rbf_full_model.mat");
% open simulink model
model = 'nonlinear_double_track_model';
% open_system(model);

%% Batch Run simulations

run_simulation = true; % warning, takes about 30 minutes

if (run_simulation)
    velocity_sweep = [1, 2, 10, 20, 30, 40]; % m/s
    steering_angle_sweep = linspace(10, 180, 5); % degress for steering wheel angle
    
    num_sims = length(velocity_sweep) * length(steering_angle_sweep);
    in(1:num_sims) = Simulink.SimulationInput(model);
    
    sim_index = 1;
    for i = 1:length(velocity_sweep)
        for u = 1:length(steering_angle_sweep)
            % run each input combination
            in(sim_index) = setBlockParameter(in(sim_index), [model '/Vehicle Plant Model/Initial Velocity (m/s)'], 'Value', num2str(velocity_sweep(i)));
            
            in(sim_index) = setBlockParameter(in(sim_index), [model '/Step Steer Input (deg)'], 'After', num2str(steering_angle_sweep(u)));
    
            sim_index = sim_index + 1;
        end
    end
    
    % Run the simulations
    out = parsim(in, 'ShowSimulationManager', 'on');

    % Save the output results to a file
    % NOTE: the data files are very large and should not be saved on the
    % Github repository. This data directory is ignored when committing
    save('./data/step_steer_sweep.mat', 'out');
else
    % Load previously saved results for analysis
    load('./data/step_steer_sweep.mat');
end

%% Plot Yaw Rate vs Time
% Number of simulations stored in out
numSims = length(out);

% Create figure for all simulation runs
figure; hold on;

% Set color order and line styles for better visibility
colors = lines(numSims); % Use default color scheme
lineStyles = {'-', '--', ':', '-.'}; % Different line styles

% Loop through each simulation run
for i = 1:numSims
    % Extract logged data (assuming logsout is used for data logging)
    logs = out(i).logsout;
    
    % Get the yaw rate signal (update the exact name)
    signal = logs.get('Yaw Rate (rad/s)'); % Replace with actual signal name
    time = signal.Values.Time;
    values = signal.Values.Data;

    % Check if any yaw rate value exceeds 3 rad/s - we can safely consider
    % these to be cars that spin out and are not steady state
    %{
    if any(values > 3)
        fprintf('Skipping Simulation %d (Yaw Rate exceeded 3 rad/s)\n', i);
        continue; % Skip this simulation
    end
    %}
    
    % Choose line style based on simulation index
    lineStyle = lineStyles{mod(i-1, length(lineStyles)) + 1};
    
    % Plot the signal with unique color and line style
    plot(time, values, 'Color', colors(i, :), 'LineWidth', 1.8, ...
        'LineStyle', lineStyle, 'DisplayName', ['Sim ' num2str(i)]);
end

% Formatting
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Yaw Rate (rad/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Yaw Rate vs Time - Base Vehicle', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Add legend to differentiate simulation runs
legend('show', 'Location', 'best', 'FontSize', 10);

% Improve appearance
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold off;