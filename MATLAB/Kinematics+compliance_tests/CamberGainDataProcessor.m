% Camber Gain Script that calculates camber gain from wheel travel and
% camber angle. It also estimates ride senstivity, optimal static camber,
% and outputs pretty graphs for design judges.

% Abigail Tucker 12/23/25

clear
clc

%replace filepath with the your file location. Use associated camber gain
%excel sheet. Remember to clsoe excel before running
data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Camber Gain\Camber Gain Excel.xlsx");

WheelTravel = data.WheelDisplacement;

CamberAngle = data.CamberAngle;


% measure distance from full droop to ride height, input value.
RideHeight = 60;

%% Camber Gain Calculation

CamberGain = gradient(CamberAngle, WheelTravel);

subplot(3, 1, 1)
plot(WheelTravel, CamberAngle,"Color", 'r', 'LineWidth',3)
xlabel('Wheel Travel [mm]')
ylabel('Camber Angle [deg]')
title('Camber Gain')
grid on
hold on

%% Camber gain vs Wheel travel
% the less linear the graph is the more senstive to setup the car will be.
% Pay attention to gain near ride height. 


subplot(3, 1, 2)
plot(WheelTravel, CamberGain, "Color",'b', 'LineWidth',3)
ax = gca;
ax.YAxis.TickLabelFormat = '%.3g';
xlabel('Wheel Travel [mm]')
ylabel('Camber Gain [deg/mm]')
grid on
hold on

yl = ylim;
% This deifines operating range based on other fsae estiates, please
% replace with better estimate if possible.
operatingRange = patch([45 75 75 45], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      [0.85 0.85 0.85], ...
      'FaceAlpha', 0.5, ...
      'EdgeColor', 'b');
xline(RideHeight, '--k', 'Ride Height', 'LineWidth', 1.5)

uistack(operatingRange,'bottom')  % moves patch behind all lines
%% Quantified RHS

operatingWindow = [45 75]; % mm
opIdx = WheelTravel >= operatingWindow(1) & WheelTravel <= operatingWindow(2);

meanGain = mean(CamberGain(opIdx));
stdGain  = std(CamberGain(opIdx)); % sensitivity metric

%% Punched in plot of operating area

subplot(3,1,3)
plot(WheelTravel(opIdx), CamberGain(opIdx))
yline(meanGain,'--k','Mean Gain')
grid on


%% Values
fprintf('\n===== CAMBER GAIN DESIGN METRICS =====\n')
fprintf('Mean Gain: %.4f\n', meanGain)
fprintf('Ride height tuning sensitivity: %.4f\n', stdGain)
fprintf('=====================================\n')






