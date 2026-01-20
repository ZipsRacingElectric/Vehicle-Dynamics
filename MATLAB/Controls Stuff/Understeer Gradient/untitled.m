%% Understeer Gradient Calculator - Clean Version

% Filepath to CSV
filepath = "C:\Users\ATuck\OneDrive - The University of Akron\Abbie's ZR Stuff\Abbie's MATLAB\Understeer Gradient Data.csv"

% Load table
Data = readtable(filepath, 'VariableNamingRule','preserve');

%% Extract steering angle and lateral acceleration
steer = Data.('CAN1.VCU_SENSOR_INPUT.STEERING_ANGLE');  % degrees
Ay = Data.('CAN1.ECUMASTER_GPS_IMU1.Y_ACCELERATION');  % m/s^2

%% Clean data - keep positive values only
% validMask = (Ay > 0) & (steer > 0);
% firstIdx = find(validMask, 1, 'first');
% T_clean = Data(firstIdx:end, :);

% steerClean = T_clean.('CAN1.VCU_SENSOR_INPUT.STEERING_ANGLE');
% AyClean = T_clean.('CAN1.ECUMASTER_GPS_IMU1.Y_ACCELERATION');

%% Smoothing
 windowSize = 10; % samples for moving average
 steerSmooth = movmean(steer, windowSize);
 AySmooth = movmean(Ay, windowSize);

%% Remove very low lateral acceleration points to avoid spikes
AyThreshold = 1; % m/s^2 (~0.1 g)
valid = AySmooth > AyThreshold;
steerSmooth = steerSmooth(valid);
AySmooth = AySmooth(valid);

%% Compute understeer gradient K
SR = 1; % Steering ratio - replace with actual
steerRad = deg2rad(steerSmooth); % convert to radians
Ay_g = AySmooth / 9.81;          % convert to g

K = (steerRad ./ Ay_g) / SR;

%% Optional: downsample for cleaner plotting
N = 5;
steerRad_ds = steerRad(1:N:end);
Ay_g_ds = Ay_g(1:N:end);
K_ds = K(1:N:end);

%% Plot results
figure;
plot(Ay_g_ds, K_ds, 'o-', 'LineWidth', 1.5)
xlabel('Lateral Acceleration [g]')
ylabel('Understeer Gradient K [rad/g]')
title('Smoothed Understeer Gradient vs Lateral Acceleration')
grid on
