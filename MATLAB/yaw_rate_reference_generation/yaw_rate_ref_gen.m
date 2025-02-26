%{

%% Overview
This script uses a quasi-static single track vehicle model to derive an ideal
yaw-rate from kinematic steering and velocity inputs. The yaw rate is valid 
for quasi-static steady state conditions, which means it is an ideal target
at corner apex, where the yaw moment is zero. This ideal yaw-rate can then 
be used as a refernce value in a yaw rate feedback control system to 
improve corner apex lateral acceleration, and thereforce overall corner
speed and lower lap times.

%% Details:
It calulates this yaw rate based on an understeer gradiant K_u which is 
designed to achieve an ideal response. Using a linear model for this is 
okay because it represents ideal vehicle behavior. By setting K_u close to
zero, the vehicle model experiences neutral steer, which gives the smallest
turning radius for a given lateral acceleration. Therefore, for a fixed
radius corner, the vehicle closest to neutral steer achieves the highest
lateral acceleration. The yaw-rate associated with this condition is
calculated and can then be the torque-vectoring controller's target.

%% TODO:
-  update cornering stiffnesses with tire data
- graph steering angle as function of a_y

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

% This understeer gradient we will specifiy to be zero, which will design
% the ideal yaw rate
K_u_ref = 0;

%% Yaw Rate vs Steering Angle
v = [10, 20, 30, 40]; % Velocities in m/s
delta_sw = linspace(0, 180, 100); % Steering wheel angle (degrees)
delta = deg2rad(delta_sw .* 0.138); % Convert steering wheel angle to dynamic wheel angle
yaw_rate = zeros(length(v), length(delta)); % yaw_rate using the linear understeer gradient
yaw_rate_ref = zeros(length(v), length(delta)); % yaw_rate using the designed understeer gradient

for i = 1:length(v)
    % Calculate maximum yaw rate
    yaw_rate_max = sigma * mu * zr25.g / v(i);

    for u = 1:length(delta)

        rate = v(i) / (zr25.wheelbase + K_u * v(i)^2) * delta(u); % Y = V / (l + K_u * V^2) * delta
        if (rate < yaw_rate_max)

            yaw_rate(i, u) = rate;
        else
            yaw_rate(i, u) = yaw_rate_max;
        end

        rate_ref = v(i) / (zr25.wheelbase + K_u_ref * v(i)^2) * delta(u); % Y = V / (l + K_u_ref * V^2) * delta
        if (rate_ref < yaw_rate_max)

            yaw_rate_ref(i, u) = rate_ref;
        else
            yaw_rate_ref(i, u) = yaw_rate_max;
        end
    end
end

%% Plot Yaw Rates
line_styles = {'-','--',':','-.'}; 
colors = lines(length(v));

figure;
hold on;
grid on;
xlabel("Steering Wheel Angle (°)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Yaw Rate (rad/s)", 'FontSize', 12, 'FontWeight', 'bold')
title("Yaw Rate vs. Steering Angle for K_u = ", K_u, ...
      'FontSize', 14, 'FontWeight', 'bold')
legend_entries = cell(1, length(v));
for i = 1:length(v)
    plot(delta_sw, yaw_rate(i, :), 'LineWidth', 2, 'LineStyle', line_styles{i}, 'Color', colors(i, :))
    legend_entries{i} = sprintf("V = %d m/s", v(i));
end
legend(legend_entries, 'Location', 'best', 'FontSize', 10)
hold off;

figure;
hold on;
grid on;
xlabel("Steering Wheel Angle (°)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Yaw Rate (rad/s)", 'FontSize', 12, 'FontWeight', 'bold')
title("Yaw Rate Reference vs. Steering Angle for K_u = ", K_u_ref, ...
      'FontSize', 14, 'FontWeight', 'bold')
legend_entries = cell(1, length(v));
for i = 1:length(v)
    plot(delta_sw, yaw_rate_ref(i, :), 'LineWidth', 2, 'LineStyle', line_styles{i}, 'Color', colors(i, :))
    legend_entries{i} = sprintf("V = %d m/s", v(i));
end
legend(legend_entries, 'Location', 'best', 'FontSize', 10)
hold off;

%% Steering Angle vs Ay

% Note: these equations to solve for steer vs Ay need to be confirmed. It
% seems like there are multiple different solutions floating around on the
% internet

% delta_dyn = delta_kin + K_u * a_y
% a_y = V^2 / R -> R = V^2 / a_y
% delta_kin = L / R = L * a_y / V^2
% delta_dyn = (L * a_y / V^2) + K_u * a_y
a_y = linspace(0, 3, 99) * zr25.g;

delta_dyn = zeros(length(v), length(a_y)); % steered angle using the linear understeer gradient
delta_dyn_ref = zeros(length(v), length(a_y)); % steered angle using the designed understeer gradient

for i = 1:length(v)
    % Calculate maximum a_y
    a_y_max = sigma * mu * zr25.g;

    for u = 1:length(a_y)

        angle = (zr25.wheelbase * a_y(u) / v(i)^2) + K_u * a_y(u);
        if (a_y(u) < a_y_max)

            delta_dyn(i, u) = (zr25.wheelbase * a_y(u) / v(i)^2) + K_u * a_y(u);
            delta_dyn_ref(i, u) = (zr25.wheelbase * a_y(u) / v(i)^2) + K_u_ref * a_y(u);
        else
            delta_dyn(i, u) = NaN;
            delta_dyn_ref(i, u) = NaN;
        end
    end
end

%% Plot Steering Angle vs Ay
line_styles = {'-','--',':','-.'}; 
colors = lines(length(v));

figure;
hold on;
grid on;
xlabel("Lateral Acceleration Ay, (m/s^2)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Steering Angle (rad)", 'FontSize', 12, 'FontWeight', 'bold')
title("Steering Angle vs A_y for K_u = ", K_u, ...
      'FontSize', 14, 'FontWeight', 'bold')
legend_entries = cell(1, length(v));
for i = 1:length(v)
    plot(a_y, delta_dyn(i, :), 'LineWidth', 2, 'LineStyle', line_styles{i}, 'Color', colors(i, :))
    legend_entries{i} = sprintf("V = %d m/s", v(i));
end
legend(legend_entries, 'Location', 'best', 'FontSize', 10)
hold off;

figure;
hold on;
grid on;
xlabel("Lateral Acceleration Ay, (m/s^2)", 'FontSize', 12, 'FontWeight', 'bold')
ylabel("Steering Angle (rad)", 'FontSize', 12, 'FontWeight', 'bold')
title("Steering Angle vs A_y for K_u_ref = ", K_u_ref, ...
      'FontSize', 14, 'FontWeight', 'bold')
legend_entries = cell(1, length(v));
for i = 1:length(v)
    plot((a_y ./ zr25.g), delta_dyn_ref(i, :), 'LineWidth', 2, 'LineStyle', line_styles{i}, 'Color', colors(i, :))
    legend_entries{i} = sprintf("V = %d m/s", v(i));
end
legend(legend_entries, 'Location', 'best', 'FontSize', 10)
hold off;