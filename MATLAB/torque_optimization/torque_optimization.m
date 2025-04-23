%{
%% Overview
Test script to validate torque optimization algorithm before implementation
into the control system

%% Details:
This takes an instantaneous snapshot of the vehicle cornering, and a yaw
moment command input. Then it solves for new motor torques that will get us
closer to the target yaw moment

%% TODO:

%}

clear; close all; clc;
Simulink.sdi.clear; % clear old simulink runs

%% Load Data
% Just using some saved data from the nonlinear_double track model
% response. V = 15 m/s, delta = 20 degrees
load("data/step_steer_sweep_single.mat");

% load QP matricies
load('QP_matrices.mat');

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

% load tire data
load("../tire_data/d2704/d2704_7in_rbf_full_model.mat");

%% Set up cornering situation
% Extract logged data (assuming logsout is used for data logging)
logs = out.logsout;
signal = get(out.logsout, 1);
vehicle_sensors = logs.get('vehicle_sensors');
slip_angles = logs.get('Slip Angles (rad)').Values;
slip_ratios = logs.get('Slip Ratios').Values;
steering_angles = logs.get('Front Steering Angles (rad)');
tire_forces = logs.get('tire_forces').Values;
yaw_ref_gen = logs.get('Yaw Rate Referece (rad/s)').Values.Data;
yaw_rate = vehicle_sensors.Values.yaw_rate.Data;
time = vehicle_sensors.Values.yaw_rate.Time;
yaw_moment_tv = logs.get('Mz_tv').Values.Data;

% Clean up extracted data into normal matricies
slip_angles = [slip_angles.Data(:,1) slip_angles.Data(:,2), slip_angles.Data(:,3), slip_angles.Data(:,4)];
slip_ratios = [slip_ratios.Data(:,1) slip_ratios.Data(:,2), slip_ratios.Data(:,3), slip_ratios.Data(:,4)];
steering_angles = [getElement(steering_angles, 2).Values.Data(:,1) getElement(steering_angles, 2).Values.Data(:,2)];
tire_data = tire_forces.Data;     % 2x4x829
numTimesteps = length(tire_forces.Time);
tire_forces = reshape(permute(tire_data, [3 1 2]), numTimesteps, 8);
tire_fy = [tire_forces(:,1) tire_forces(:,3) tire_forces(:,5) tire_forces(:,7)];
tire_fx = [tire_forces(:,2) tire_forces(:,4) tire_forces(:,6) tire_forces(:,8)];

% Plot transient response of run
figure;
plot(time, yaw_ref_gen, 'r--', 'LineWidth', 1.5); % Yaw reference (dashed red line)
hold on;
plot(time, yaw_rate, 'b-', 'LineWidth', 1.5); % Yaw rate (solid blue line)
xlabel('Time (s)');
ylabel('Yaw (deg/s)');
title('Yaw Reference vs. Yaw Rate');
legend('Yaw Reference', 'Yaw Rate');
grid on;
hold off;

%% Grab instanteous data for testing torque optimization
index = 55;
t = time(index);
apps = 0;

% 1st objective we will try to achieve
% We will take the derivative of the measured yaw rate sensor
time_step = time(index) - time(index - 1);
[yaw_moment, ~] = get_yaw_moment(yaw_rate(index), yaw_rate(index - 1), time_step, zr25.yaw_polar_inertia);
Mz_tv = yaw_moment_tv(index);

% 2nd object we will try to achieve
x_dot = vehicle_sensors.Values.body_acceleration_x.Data;
ax_ref = get_x_accel_target(apps);

%% Calculate VD Objectives and other values

% Slip Angles
alpha_1 = slip_angles(index, 1);
alpha_2 = slip_angles(index, 2);
alpha_3 = slip_angles(index, 3);
alpha_4 = slip_angles(index, 4);

% Slip Ratios
s_fl_0 = slip_ratios(index, 1);
s_fr_0 = slip_ratios(index, 2);
s_rl_0 = slip_ratios(index, 3);
s_rr_0 = slip_ratios(index, 4);

% Steering Angles
delta_fl = steering_angles(index, 1);
delta_fr = steering_angles(index, 2);
delta_rl = 0;
delta_rr = 0;

% Tire Forces
Fx_fl_0 = tire_fx(index, 1);
Fx_fr_0 = tire_fx(index, 2);
Fx_rl_0 = tire_fx(index, 3);
Fx_rr_0 = tire_fx(index, 4);

Fy_fl_0 = tire_fy(index, 1);
Fy_fr_0 = tire_fy(index, 2);
Fy_rl_0 = tire_fy(index, 3);
Fy_rr_0 = tire_fy(index, 4);

% Tire Force derivatives

fy = get_sl_slice(Fy_grid, sl_vals, sa_vals, fz_vals, -4, -600);
figure;
plot(sl_vals, fy);
xlabel('Slip Ratio');
ylabel('Interpolated Fy');
title('Interpolated Fy vs Slip Ratio at sa=3.5°, Fz=2750 N');

fx = get_sl_slice(Fx_grid, sl_vals, sa_vals, fz_vals, -4, -600);
figure; 
plot(sl_vals, fx);
xlabel('Slip Ratio');
ylabel('Interpolated Fx');
title('Interpolated Fx vs Slip Ratio at sa=3.5°, Fz=2750 N');

dFy_fl_ds = get_derivative(sl_vals, fy, s_fl_0);
dFy_fr_ds = get_derivative(sl_vals, fy, s_fr_0);
dFy_rl_ds = get_derivative(sl_vals, fy, s_rl_0);
dFy_rr_ds = get_derivative(sl_vals, fy, s_rr_0);

dFx_fl_ds = get_derivative(sl_vals, fx, s_fl_0);
dFx_fr_ds = get_derivative(sl_vals, fx, s_fr_0);
dFx_rl_ds = get_derivative(sl_vals, fx, s_rl_0);
dFx_rr_ds = get_derivative(sl_vals, fx, s_rr_0);

%{
plot_linear(sl_vals, fy, dFy_fl_ds, sl_fl_0);
plot_linear(sl_vals, fy, dFy_fr_ds, sl_fr_0);
plot_linear(sl_vals, fy, dFy_rl_ds, sl_rl_0);
plot_linear(sl_vals, fy, dFy_rr_ds, sl_rr_0);

plot_linear(sl_vals, fx, dFx_fl_ds, sl_fl_0);
plot_linear(sl_vals, fx, dFx_fr_ds, sl_fr_0);
plot_linear(sl_vals, fx, dFx_rl_ds, sl_rl_0);
plot_linear(sl_vals, fx, dFx_rr_ds, sl_rr_0);
%}


%% Calculate QP Matrix inputs
% The optimization needs to be constrained around these limits

% weights and scaling factors
w1 = 0.60;
w2 = 0.25;
w3 = 0.10;
tf = zr25.track_width_front;
tr = zr25.track_width_rear;
mass = zr25.mass_total;
lf = zr25.a;
lr = zr25.b;
R = zr25.tire_loaded_radius;

%% Solve QP Matrix
params.A = [1; 1; 1; 1]; % Constraint matrix
params.P = double(simplify(subs(P))); % Quadratic matrix
params.l = [1; 1; 1; 1]; % Matrix of minimum constraints
params.q = double(simplify(subs(q))); % Linear matrix
params.u = [1; 1; 1; 1]; % Matrix of maximum constraints

%% Run Optimization

%}
%% Helper Functions
% Calculate yaw moment from derivative of yaw rate sensor
function [yaw_moment, yaw_acceleration]  = get_yaw_moment(yaw_rate, yaw_rate_prior, time_step, inertia)
    yaw_acceleration = yaw_rate - yaw_rate_prior / time_step; % rad/s^2
    yaw_moment = yaw_acceleration * inertia;
end

% Calculate longitudinal acceleration target from driver sensors
function ax_ref  = get_x_accel_target(apps)
    % placeholder for some more complicated algorithm in the future
    ax_max = 2 * 9.81; % m/s^2
    ax_ref =  ax_max * apps / 100;
end

% Grab a tire_force vs slip ratio at a slip angle and normal force operating point
function sl_slice = get_sl_slice(tire_grid, sl_vals, sa_vals, fz_vals, sa_query, fz_query)
    % Preallocate output slice
    sl_slice = zeros(size(sl_vals));
    
    % Loop over slip ratio values and interpolate for each one
    for i = 1:length(sl_vals)
        sl = sl_vals(i);
        sl_slice(i) = interpn( ...
            sl_vals, sa_vals, fz_vals, tire_grid, ...
            sl, sa_query, fz_query, 'linear');
    end
end

% Calculate the deriviative at a given operating point
function dF_ds = get_derivative(sl_vals, Fy_vals, sl_point)
    % Use spline interpolation to get smooth derivative
    pp = spline(sl_vals, Fy_vals);           % Piecewise polynomial
    pp_deriv = fnder(pp);                    % Derivative of the spline
    dF_ds = ppval(pp_deriv, sl_point);     % Evaluate derivative at query point
end

% Plots the tire data along with its linearized version
function plot_linear(sl_vals, Fy_vals, dF_ds, sl_query)
    % Get interpolated Fy and derivative at query point
    Fy_query = interp1(sl_vals, Fy_vals, sl_query, 'spline');
    slope = dF_ds;
    
    % Create tangent line around the query point
    sl_range = linspace(min(sl_vals), max(sl_vals), 200);
    tangent_line = Fy_query + slope * (sl_range - sl_query);

    % Plot everything
    figure;
    plot(sl_vals, Fy_vals, 'b-', 'LineWidth', 2); hold on;
    plot(sl_query, Fy_query, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(sl_range, tangent_line, 'r--', 'LineWidth', 1.5);
    xlabel('Slip Ratio');
    ylabel('Lateral Force Fy');
    legend('Fy vs Slip', 'Operating Point', 'Tangent Line');
    title(sprintf('Tangent at Slip = %.3f', sl_query));
    grid on;
end