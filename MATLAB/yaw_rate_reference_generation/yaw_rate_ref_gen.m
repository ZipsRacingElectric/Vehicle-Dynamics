%{

%% Overview
This script uses a linear single track vehicle model to derive an ideal
yaw-rate from kinematic steering and velocity inputs. This ideal yaw-rate
can then be used as a refernce value in a yaw rate feedback control system.

%% Details:
It calulates this yaw rate based on an understeer gradiant K_u which is designed to achieve an
ideal response. By setting K_u close to zero, the vehicle model experiences
neutral steer, which gives the smallest turning radius for a given lateral
acceleration. Therefore, for a fixed radius corner, the vehicle closest to
neutral steer achieves the highest lateral acceleration. The yaw-rate
associated with this condition is calculated and can then be the
torque-vectoring controler's target.

%% TODO:
-  update cornering stiffnesses with tire data
- graph yaw rate using the target understeer gradient

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

%% Define linear vehicle understeer gradient
C_f = 0.18; % front cornering stiffness placeholder
C_r = 0.17; % front cornering stiffness placeholder
Fz_f = zr25.front_mass * zr25.g;
Fz_r = zr25.rear_mass * zr25.g;
sigma = 2/3; % tire-road friction modifier
mu = 3; % max tire-road friction coefficient from TTC data

% As expected, this is slightly oversteering for a linear FSAE vehicle
% without compliance
K_u = (zr25.mass_total / zr25.wheelbase) * ((zr25.b / (C_f * Fz_f)) -  (zr25.a / (C_r * Fz_r)));

%% Graph Yaw Rate
v = [10, 20, 30, 40]; % Velocities in m/s
delta_sw = linspace(0, 180, 100); % Steering wheel angle (degrees)
delta = deg2rad(delta_sw .* 0.138); % Convert steering wheel angle to dynamic wheel angle

line_styles = {'-','--',':','-.'}; 
colors = lines(length(v));
figure;
hold on;
grid on;
xlabel("Steering Wheel Angle (°)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Yaw Rate (rad/s)", 'FontSize', 12, 'FontWeight', 'bold')
title("Yaw Rate vs. Steering Wheel Angle for Linear Understeer Gradient", ...
      'FontSize', 14, 'FontWeight', 'bold')
legend_entries = cell(1, length(v));

for i = 1:length(v)
    
    % Calculate maximum yaw rate
    yaw_rate_max = sigma * mu * zr25.g / v(i);

    % Calculate yaw rate
    yaw_rate = zeros(1, length(delta));
    for u = 1:length(delta)
        rate = delta(u) / (zr25.wheelbase + K_u * v(i)); % Y = delta_dyn / (l + K_u * V)
        if (rate < yaw_rate_max)
            yaw_rate(u) = rate;
        else
            yaw_rate(u) = yaw_rate_max;
        end
    end


    plot(delta_sw, yaw_rate, 'LineWidth', 2, 'LineStyle', line_styles{i}, 'Color', colors(i, :))
    legend_entries{i} = sprintf("V = %d m/s", v(i));
end

legend(legend_entries, 'Location', 'best', 'FontSize', 10)
hold off;
