% Bump Steer Test Data Process
% This code takes data from Bump Steer Testing and calculate bump steer
% angle
%Written by Carl and Khoi
Excel = readtable("/Users/khoi/Library/CloudStorage/OneDrive-TheUniversityofAkron/200 Controls/Data Excels/bumpSteer_Eliminator.xlsx")

% Importing all data

Z = Excel.ZMovement(1:10) % Data for Independent Z movement

D = Excel.distanceOnWall(1:10) % Data for Dependent D change

dial = Excel.dialIndicator(1:10) % Data for Dependent dial indicator

Y = Excel.distanceToWall(1:10) % Data for Fixed Y


% Calculate X using distanceOnWall and dialIndicator

X = D - dial


% Calculate bump steer angle

bumpSteerAngle = atand(X./Y) % Calculate bump steer angle in degrees


% Visual Graph for bumpSteerAngle

x = Z
y = bumpSteerAngle
plot(x, y)
xlabel('Z Movement')
ylabel('Bump Steer Angle')
title('Bump Steer Angle vs Z Movement');
grid on;
