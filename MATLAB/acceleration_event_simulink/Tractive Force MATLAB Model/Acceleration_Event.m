% Parameters
mass = 300; % mass of the car in kg
motorTorque = 1000; % torque of the electric motor in Nm
finalDriveRatio = 3.5; % final drive ratio of the car
tireRadius = 0.25; % loaded radius of the tire in meters
distance = 90; % distance of the acceleration run in meters

% Calculate acceleration
time = sqrt((2 * distance) / (motorTorque * finalDriveRatio * tireRadius * mass));

% Create time vector
t = linspace(0, time, 1000);

% Calculate tractive force
Ft = motorTorque * finalDriveRatio / tireRadius;


% Plot tractive force over time
plot(t, Ft * ones(size(t)), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Tractive Force (N)');
title('Tractive Force over Time');
grid on;
hold on;

% Plot acceleration
plot(t, Ft * ones(size(t)), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Tractive Force (N)');

% Plot acceleration

%Plot Power of motor


