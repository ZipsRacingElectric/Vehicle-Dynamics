

%% Main Jawn
% Calculates camber and toe angle change for defined forces and moments
% still needs data for spring constants of knuckles and rod ends
% This will be used to compare real test data from the car for design
% purposes

% Written by Abbie Tucker, spring 2025


% Define file paths
filepath = 'C:\Users\benmo\OneDrive - The University of Akron\Documents - Zips Racing FSAE\ZR25\Vehicle Dynamics\System\Analysis\'
pointsfile = 'ZR25_SuspensionForcesHeavy.xlsx'

% Load suspension points
FL = SuspensionPoints_better(fullfile(filepath, pointsfile), false);
disp(class(FL))

% Define suspension point names
FLpointPairs =getPointPairs(FL);

% Uses position and Unit Vectors to assemble a matrix containing a system
% of 6 equations.
SystemMatrix = BuildMatrix(FL, FLpointPairs);

% Suspension Component Data
% ordered as:R_LCA_FORE, R_LCA_AFT, R_UCA_FORE, R_UCA_AFT, R_TIE, R_PUSH
GetSuspensionData

% Define Forces and Moments
TestCase = [0; 0; 0.0; 5.3761; 6.5259; 0] * 1e5;

% Calculates spring constants
Ksus = SpringySpring(GetSuspensionData);
%Kre = SpringySpring(GetSuspensionData);
Ks = [Ksus];
SusSpringSeries = SpringSeries(Ks);

% Calculates and Orders Forces
Forces = solvySolve(SystemMatrix, TestCase);
ReorderedForces  = [Forces.LCAFore Forces.LCAAft Forces.UCAFore Forces.UCAAft Forces.Toe Forces.Push];

% Calculates Displacement and Displacement Vectors
Displacements = getDisplacements(SusSpringSeries, ReorderedForces, getUnitVectors(FL,FLpointPairs));

% Calculate Camber Angle
CamberAnglePoints = GetTotalCDisplacement(Displacements);
ZDistance = getZDistance(FL.ucaOutboard(3,1), FL.lcaOutboard(3,1));
CamberAngleChange = getAngle(CamberAnglePoints,ZDistance);


% Calculate Toe Angle 
ToeAnglePoints = GetTotalTDisplacement(Displacements);
XDistance = getXDistance(FL.ucaOutboard(1,1),FL.toeOutboard(1,1));
ToeAngleChange = getAngle(ToeAnglePoints,XDistance);

%% Plotting Camber and Toe Angle Changes vs Applied Moment

% Define a range of moments to apply (for example, varying only Fx or Mz)
momentSweep = 0:10000:30000;  % Apply Mx from 0 to 150 N·m in 25 N·m increments
numTests = length(momentSweep);  % Update the number of test cases

camberChanges = zeros(1, numTests);
toeChanges = zeros(1, numTests);

for i = 1:numTests
    TestCase = [0; 0; 0.0089; momentSweep(i); 0; 0]; % Only Mx is changing
    Forces = solvySolve(SystemMatrix, TestCase);
    ReorderedForces  = [Forces.LCAFore Forces.LCAAft Forces.UCAFore Forces.UCAAft Forces.Toe Forces.Push];
    Displacements = getDisplacements(SusSpringSeries, ReorderedForces, getUnitVectors(FL,FLpointPairs));
    
    CamberAnglePoints = GetTotalCDisplacement(Displacements);
    ZDistance = getZDistance(FL.ucaOutboard(3,1), FL.lcaOutboard(3,1));
    camberChanges(i) = getAngle(CamberAnglePoints, ZDistance);

    ToeAnglePoints = GetTotalTDisplacement(Displacements);
    XDistance = getXDistance(FL.ucaOutboard(1,1),FL.toeOutboard(1,1));
    toeChanges(i) = getAngle(ToeAnglePoints, XDistance);
end

% Plotting
figure;
%plot(momentSweep / 1000, camberChanges, '-s', 'LineWidth', 2); 
hold on;
plot(momentSweep/1000, toeChanges*100);
xlabel('Applied Mx (N·m)');
ylabel('Deflection (degrees)');
title('Toe Deflection vs Applied Mx');
legend('Predicted Toe Change', 'Measured Toe Change');
grid on;
set(gca, 'FontSize', 12);

knownMoments = .001*9.81.*[0, 2000, 3000, 4000]*.5;  % N·mm
knownDisplacements = rad2deg(atan([0, 4, 6.5, 9]./1473)); % mm displacement 

% Plot

plot(knownMoments, knownDisplacements, 'o-', 'LineWidth', 2);
% xlabel('Applied Moment (N·m)');
% ylabel('Measured Toe Displacement (mm)');
% title('Measured Displacement vs Applied Moment');
% grid on;
% set(gca, 'FontSize', 12);
legend('Predicted Toe Change', 'Measured Toe Change');

%% Functions

function PointPairs = getPointPairs(~)

 PointPairs = {
        'pushInboard', 'pushOutboard';
        'ucaFore', 'ucaOutboard';
        'ucaAft', 'ucaOutboard';
        'lcaFore', 'lcaOutboard';
        'lcaAft', 'lcaOutboard';
        'toeInboard', 'toeOutboard'}
end

function SuspensionDataStruct = GetSuspensionData()
    
   SuspensionDataStruct.lengths = [287.05 277.27 255.42 271.54 275.28 306.55]
   SuspensionDataStruct.tubeODs = [0.625 0.625 0.5 0.5 0.5 0.5].*25.4
   SuspensionDataStruct.tubeWallThiccs = [0.028 0.028 0.028 0.035 0.028 0.028].*25.4
   SuspensionDataStruct.youngsModi = 205000*ones(size(SuspensionDataStruct.lengths));

end

function tubeArea = getTubeArea(tubeOD,tubeWallThicc)
    tubeArea = .25*pi*(tubeOD.^2 - (tubeOD - 2.*tubeWallThicc ).^2);
end

function UnitVectors = getUnitVectors(suspensionStuff, pointPairs)
    
    % Compute unit vectors
    numPoints = size(pointPairs, 1);
    UnitVectors = zeros(3, numPoints);
  
    
    for i = 1:numPoints

        inboard = suspensionStuff.(pointPairs{i, 1});
        outboard = suspensionStuff.(pointPairs{i, 2});
        
        UnitVectors(:, i) = (outboard - inboard) / norm(outboard - inboard);
    end
end

function PositionMatrix = getPositionMatrix(suspensionStuff, pointPairs)
     numPoints = size(pointPairs, 1);
    PositionMatrix = zeros(3, numPoints);
    
    for i = 1:numPoints
        inboard = suspensionStuff.(pointPairs{i, 1});
        outboard = suspensionStuff.(pointPairs{i, 2});

        PositionMatrix(:, i) = inboard; % Store inboard points for cross-product matrix
    end
end

function SystemMatrix = BuildMatrix(suspensionStuff, pointPairs)

    UnitVectors = getUnitVectors(suspensionStuff, pointPairs)
    PositionMatrix = getPositionMatrix(suspensionStuff, pointPairs)

   % Construct matrices
    UxPMatrix = cross(PositionMatrix, UnitVectors);
    SystemMatrix = [UnitVectors; UxPMatrix];
end

function ForcesStruct = solvySolve(SystemMatrix,knownForceMomentThing)    
    
    %% the solvysolve part
    Forces = linsolve(SystemMatrix, knownForceMomentThing);
    
    
    
    ForcesStruct.Push = Forces(1,1)
    ForcesStruct.UCAFore = Forces(2,1)
    ForcesStruct.UCAAft = Forces(3,1)
    ForcesStruct.LCAFore = Forces(4,1)
    ForcesStruct.LCAAft = Forces(5,1)
    ForcesStruct.Toe = Forces(6,1)

end

function SpringConstant = SpringySpring(SuspensionDataStruct)

    SpringConstant = (SuspensionDataStruct.youngsModi.*getTubeArea(SuspensionDataStruct.tubeODs, SuspensionDataStruct.tubeWallThiccs))./SuspensionDataStruct.lengths;

end

function SpringsInSeries = SpringSeries(Ks)
     % k's = vector of spring constants
    [rows, cols] = size(Ks); % Get the matrix dimensions

   if rows == 1 % check if Ks is (1 x n)
      SpringsInSeries = Ks

   else
      SpringsInSeries = 1 ./ sum(1./Ks)
   end

end

function Displacements = getDisplacements(Ks,Forces,unitVectors)

Displacements = Forces./Ks

   Displacements = Displacements.*unitVectors

end

function ZDistance =  getZDistance(ZpointUCA, ZPointLCA)

ZDistance = ZpointUCA-ZPointLCA

end

function XDistance = getXDistance(XpointUCA, XPointToe)

XDistance = XpointUCA-XPointToe

end

function TotalCamberDisplacement = GetTotalCDisplacement(Displacement)

    LCAADisplacementY = Displacement(2,2)
    LCAFDisplacementY = Displacement(2,1)
    UCAADisplacementY = Displacement(2,4)
    UCAFDisplacementY = Displacement(2,3)

    LCADisplacementY = (LCAADisplacementY+LCAFDisplacementY)
    UCADisplacementY = (UCAADisplacementY+UCAFDisplacementY)

    TotalCamberDisplacement = (UCADisplacementY-LCADisplacementY)
end

function TotalToeDisplacement = GetTotalTDisplacement(Displacement)


ToeDisplacement = Displacement(2,5)
 UCAADisplacement = Displacement(2,4)
 UCAFDisplacement = Displacement(2,3)
 UCADisplacement = (UCAADisplacement-UCAFDisplacement)

 TotalToeDisplacement = UCADisplacement-ToeDisplacement

end

function DegreeAngleCalculator = getAngle(DisplacementPoints, ComponentDistance)

DegreeAngleCalculator = atan(DisplacementPoints/ComponentDistance)* (180 / pi)

end
