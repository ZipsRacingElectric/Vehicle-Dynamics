
% ------------------------------------------------------------
% This MATLAB script sets up a Simulink model architecture for your
% FSAE vehicle, using parameters from your 'zr' object.
% It creates the high-level signal flow and subsystem placeholders
% for you to fill in with equations and models later.
% Abigail Tucker 10/27/25
% ------------------------------------------------------------


addpath vehicle_data ;
githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr26_data.xlsx');
Zr26 = vehicle(parameterSpreadsheet);

%zr26 = vehicle(parameterSpreadsheet);
%% Ensure zr object is available
if ~exist('Zr26', 'var')
    % Path to your meme image
    memePath = fullfile(pwd, 'SillyGoose.png'); % or .png, etc.

    if isfile(memePath)
        % Read and display the meme
        img = imread(memePath);
        figure('Name','Silly Goose Detected','NumberTitle','off');
        imshow(img);
        title('Silly Goose 🪿 – Vehicle object not found!');
    else
        warning('No meme image found at %s', memePath);
    end

    % Print the honk message
    fprintf('\n🪿 HONK! Silly goose – vehicle object ''zr'' not found!\n');
    fprintf('Hint: run vehicle.m to define it.\n\n');

    % Throw the actual MATLAB error so execution stops
    error('Silly goose 🪿, Vehicle object ''zr'' not found in workspace. Please define it first. hint: use vehicle.m');
end

%% Create Simulink model
modelName = 'ZR_ControlsArchitecture';
if bdIsLoaded(modelName)
    close_system(modelName,0);
end
new_system(modelName);
open_system(modelName);

%% Define Subsystems
subsystems = {
    'DriverCmd'               [100 100 250 150];
    %'TorqueVectorController'  [300 100 470 150];
    'PowertrainMotor'         [300 100 470 150];
    'TireModel'               [520 100 650 150];
    'VehicleDynamics'         [700 100 850 150];
    %'SensorsIMU'              [1100 100 1220 150];
};

for i = 1:size(subsystems,1)
    add_block('simulink/Ports & Subsystems/Subsystem', [modelName '/' subsystems{i,1}], 'Position', subsystems{i,2});
end

%% Add Inputs & Outputs
add_block('simulink/Sources/In1', [modelName '/SteeringAngle δ'], 'Position', [30 100 60 120]);
add_block('simulink/Sources/In1', [modelName '/PedalPos'], 'Position', [30 140 60 160]);
add_block('simulink/Sinks/Out1', [modelName '/Vx'], 'Position', [1250 100 1280 120]);
add_block('simulink/Sinks/Out1', [modelName '/Vy'], 'Position', [1250 130 1280 150]);
add_block('simulink/Sinks/Out1', [modelName '/YawRate'], 'Position', [1250 160 1280 180]);

%% Connect Lines
add_line(modelName, 'PedalPos/1', 'DriverCmd/1');
add_line(modelName, 'SteeringAngle δ/2', 'DriverCmd/2');
add_line(modelName, 'DriverCmd/1', 'PowertrainMotor/1');
add_line(modelName, 'DriverCmd/2', 'PowertrainMotor/2');
add_line(modelName, 'PowertrainMotor/1', 'TireModel/1');
add_line(modelName, 'TireModel/1', 'VehicleDynamics/1');
add_line(modelName, 'VehicleDynamics/1', 'SensorsIMU/1');
add_line(modelName, 'SensorsIMU/1', 'Vx/1');
add_line(modelName, 'SensorsIMU/2', 'Vy/1');
add_line(modelName, 'SensorsIMU/3', 'YawRate/1');

%% Annotate model
annotationText = sprintf(['ZR Vehicle Model Template\n', ...
    'Created: %s\n', ...
    'Subsystems: DriverCmd, TVC, Powertrain, TireModel, VehicleDynamics, SensorsIMU'], datestr(now));
add_block('simulink/Notes/Note', [modelName '/ModelInfo'], ...
    'Position', [100 200 400 300], 'Text', annotationText);

%% Save system
save_system(modelName);

fprintf('✅ Simulink model "%s" created successfully.\n', modelName);
fprintf('You can now populate each subsystem with equations using zr parameters.\n');