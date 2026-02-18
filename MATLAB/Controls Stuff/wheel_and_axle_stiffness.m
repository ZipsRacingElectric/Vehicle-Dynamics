clear, clc

run('/Users/khoi/Documents/GitHub/Vehicle-Dynamics/MATLAB/Kinematics+compliance_tests/MotionRatioDataProcessor.m');
%Given Motion Ratio, Spring Rate, ARB

ZR25 = vehicle('/Users/khoi/Documents/GitHub/Vehicle-Dynamics/MATLAB/vehicle_data/vehicle.m')

frontTrackWidth = 20
rearTrackWidth = 22

frontARBStiffness = 2
rearARBStiffness = 2

springRateF = [0.4, 0.5, 0.1, 10, 40, 25, 70, 150, 500, 1200];
springRateR = [0.4, 0.5, 0.1, 10, 40, 25, 70, 150, 500, 1200];

AMRF = (AMRFR + AMRFL) ./ 2
AMRR = (AMRRR + AMRRL) ./ 2


wheelRateF = springRateF .* AMRF.^2;
wheelRateR = springRateR .* AMRR.^2;

axleStiffnessF = wheelRateF .* frontTrackWidth ./ 2;
axleStiffnessR = wheelRateR .* frontTrackWidth ./ 2;

axleStiffnessF = axleStiffnessF + frontARBStiffness;
axleStiffnessR = axleStiffnessR + rearARBStiffness;

% Calculate total stiffness for the front and rear axles

A = numel(springRateF)
B = numel(frontARBStiffness)
table = nan(A,B);

% Populate the table with calculated stiffness values
for i = 1:A
    for j = 1:B
      table(i,j) = totalStiffness(axleStiffnessF,axleStiffnessR)
    end

end
     

function [avg] = totalStiffness(axleStiffnessF, axleStiffnessR) ;

result = axleStiffnessF + axleStiffnessR;
avg = mean(result);

end