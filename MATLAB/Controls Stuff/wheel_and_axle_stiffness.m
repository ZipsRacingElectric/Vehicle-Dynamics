clear, clc
%Given Wheel Rate, Spring Rate, ARB
addpath vehicle_data ;
githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr25_data.xlsx');
ZR26 = vehicle(parameterSpreadsheet);

frontTrackWidth = ZR26.track_width_front;
rearTrackWidth = ZR26.track_width_rear;

rearARBStiffness = ZR26.wheel_rate_from_bar_front;
springRate = [0.4, 0.5, 0.1, 10, 40, 25, 70, 150, 500, 1200];

wheelRateF = ZR26.wheel_rate_front;
wheelRateR = ZR26.wheel_rate_rear;

axleStiffnessF = wheelRateF .* frontTrackWidth ./ 2;
axleStiffnessR = wheelRateR .* frontTrackWidth ./ 2;

axleStiffnessR = axleStiffnessR + rearARBStiffness;


% Calculate total stiffness for the front and rear axles
rowForTable = numel(springRate);
preTableF = [springRate', repmat(wheelRateF', rowForTable, 1)];
preTableR = [springRate', repmat(rearARBStiffness', rowForTable, 1), repmat(wheelRateR', rowForTable, 1)];
axleStiffnessFTable = nan(rowForTable,1);
axleStiffnessRTable = nan(rowForTable,1);
totalRollStiffnessTable = nan(rowForTable,1);
% Populate the table with calculated stiffness values
for i = 1:rowForTable
axleStiffnessFTable(i,1) = wheelRateF .* frontTrackWidth ./ 2;
end
for i = 1:rowForTable
axleStiffnessRTable(i,1) = wheelRateR .* rearTrackWidth ./ 2 + rearARBStiffness;
end
for i = 1:rowForTable
totalRollStiffnessTable(i,1) = axleStiffnessFTable(i) + axleStiffnessRTable(i);
end

tableF = [preTableF, axleStiffnessFTable, totalRollStiffnessTable]
tableR = [preTableR, axleStiffnessRTable, totalRollStiffnessTable]
