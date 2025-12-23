% Cameber Gain Script that calculates camber gain from wheel travevl and
% camber angle. It also estimates ride senstivity, optimal static camber,
% and outputs pretty graphs for design judges.

% Abigail Tucker 12/23/25


%replace filepath with the your file location. Use associated camber gain
%excel sheet. Remember to clsoe excel before running
data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Camber Gain\Camber Gain Excel.xlsx");

WheelTravel = data.WheelDisplacement;

CamberAngle = data.CamberAngle;


% measure distance from full droop to ride height, input value.
RideHeight = 60

%% Camber Gain Calculation

CamberGain = gradient(CamberAngle, WheelTravel);

subplot(3, 1, 1)
plot(WheelTravel, CamberAngle,"Color", 'r', 'LineWidth',3)
xlabel('Wheel Travel')
ylabel('Camber Angle')
title('Camber Gain')
grid on
hold on

%% Camber gain vs Wheel travel
% the less linear the graph is the more senstive to setup the car will be.
% Pay attention to gain near ride height. 


subplot(3, 1, 2)
plot(WheelTravel, CamberGain, "Color",'b', 'LineWidth',3)
xlabel('Wheel Travel')
ylabel('Camber Gain')
grid on
hold on

yl = ylim
% This deifines operating range based on other fsae estiates, please
% replace with better estimate if possible.
operatingRange = patch([45 75 75 45], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.85 0.85 0.85], ...
      'FaceAlpha', 0.5, ...
      'EdgeColor', 'b');
xline(RideHeight, '--k', 'Ride Height', 'LineWidth', 1.5)

uistack(operatingRange,'bottom')  % moves patch behind all lines



