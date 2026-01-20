% Emily Berdis 1/16/2026
clear
clc
close all

%replace filepath with the your file location. Use associated steer angles
%excel sheet. Remember to close excel before running
% Load the data from the specified Excel file
data = readtable("C:\Users\etber\Documents\GitHub\Vehicle-Dynamics\MATLAB\Kinematics+compliance_tests\Steer-Steer Test.xlsx"); % Specify the correct file path
RWA = data.RightWheelAngle;
LWA = data.LeftWheelAngle;
SWA = data.SteeringWheelAngle;
%Plot variables
figure;
plot(RWA, 'r', 'DisplayName', 'Right Wheel Angle');
hold on;
plot(LWA, 'b', 'DisplayName', 'Left Wheel Angle');
plot(SWA, 'g', 'DisplayName', 'Steering Wheel Angle');
legend show;
xlabel('Test Count');
ylabel('Wheel Angle (degrees)');
title('Wheel Angles vs Test Count');
xticks([0,1,2,3,4,5])
grid on;

% Plot the Left and Right Wheel Angles compared to the Steering Wheel Angle
figure;
plot(SWA, RWA, 'ro-', 'DisplayName', 'Right Wheel Angle vs Steering Angle');
hold on;
plot(SWA, LWA, 'bo-', 'DisplayName', 'Left Wheel Angle vs Steering Angle');
legend show;
xlabel('Steering Wheel Angle (degrees)');
ylabel('Wheel Angle (degrees)');
title('Wheel Angles Compared to Steering Wheel Angle');
grid on;
