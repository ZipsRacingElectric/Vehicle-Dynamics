% Emily Berdis 1/16/2026
clear
clc
close all

%replace filepath with the your file location. Use associated steer angles
%excel sheet. Remember to close excel before running
% Load the data from the specified Excel file
data = readtable("C:\Users\etber\OneDrive - The University of Akron\Zips Racing FSAE - Documents\ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Steer Steer\Steer-Steer Test.xlsx"); % Specify the correct file path
RWA = data.RightWheelAngle;
LWA = data.LeftWheelAngle;
SWA = data.SteeringWheelAngle;
valid= LWA~=0; %logical index of non-zero elements
valid= RWA~=0; %logical index of non-zero elements
% Calculate the Steering Ratio Left
SRL= mean(SWA(valid) ./LWA(valid));
% Calculate the Steering Ratio Right
SRR = mean(SWA(valid) ./ RWA (valid)); % Element-wise division to compute steering ratio

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
txt= ['Steering Ratio Left =' num2str(SRL, '%.2f') ]
text(0, -15, txt, 'FontSize', 16)
txt= ['Steering Ratio Right =' num2str(SRR, '%.2f') ]
text(0, -20, txt, 'FontSize', 16)
legend show;
xlabel('Steering Wheel Angle (degrees)');
ylabel('Wheel Angle (degrees)');
title('Wheel Angles Compared to Steering Wheel Angle');
grid on;
