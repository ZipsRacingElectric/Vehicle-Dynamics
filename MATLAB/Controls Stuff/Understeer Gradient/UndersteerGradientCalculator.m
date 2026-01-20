%% Understeer Gradient Calculator
%The understeer gradient is a key metric in vehicle dynamics that describes how a car’s steering response changes with lateral acceleration. It quantifies whether a vehicle tends to understeer (requires more steering angle as cornering force increases), oversteer (requires less steering angle), or is neutral steer.
%A positive understeer gradient means the car needs progressively more steering input to hold higher lateral acceleration → typical for production cars (stable, predictable).
%A negative gradient indicates oversteer → the car needs less steering input as cornering builds, which can become unstable at the limit.
%Zero gradient is neutral steer → steering angle vs. lateral acceleration is linear.
%This metric is useful because it:
% Provides a simple, quantitative way to describe handling balance.
%Helps engineers compare setup changes (springs, bars, tire pressures, aero balance).
%Links test data (steering, lateral accel, speed) to fundamental handling theory.
%Serves as a tuning target — FSAE teams often aim for a mild understeer gradient for stability while still keeping sharp turn-in.

%driving equatioins: For constant radius/Variable Speed - K =
%(SteerAngle/LatAccel)

%Constant speed (Throttle)/ Variable Radius -
%(SteerAngle/LatAccel)-(Wheelbase/velocity^2)

% Filepath to CSV
filepath = "C:\Users\ATuck\OneDrive - The University of Akron\Abbie's ZR Stuff\Abbie's MATLAB\Understeer Gradient Data.csv"

% Load table
Data = readtable(filepath, 'VariableNamingRule','preserve');

%% Extract steering angle and lateral acceleration
steer = Data.('CAN1.VCU_SENSOR_INPUT.STEERING_ANGLE');  % degrees
Ay = Data.('CAN1.ECUMASTER_GPS_IMU1.Y_ACCELERATION');  % m/s^2

K = (steerRad ./ Ay_g) / SR;

%% Plot results
figure;
plot(Ay_g_ds, K_ds, 'o-', 'LineWidth', 1.5)
xlabel('Lateral Acceleration [g]')
ylabel('Understeer Gradient K [rad/g]')
title('Smoothed Understeer Gradient vs Lateral Acceleration')
grid on
