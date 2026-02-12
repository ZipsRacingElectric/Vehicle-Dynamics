%% Motion Ratio 

% this code will plot vertical wheel travel against corresponding damper
% length providing the motion ratio for each corner 

% Written by Abigail Tucker 01/09/26

clc, clear
close all


filePath = "C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Motion Ratio\Motion Ratio Excel.xlsx";
excel = readtable(filePath);

% Load corner Data
FrontRight = readtable(filePath, "Range", "FrontRight");
FrontLeft = readtable(filePath, "Range", "FrontLeft");
RearRight = readtable(filePath, "Range", "RearRight");
RearLeft = readtable(filePath, "Range", "RearLeft");

FRWheelDisplacment = FrontRight.WheelDisplacement_mm_;
FRDamperLength = FrontRight.DamperLength_mm_;

FLWheelDisplacment = FrontLeft.WheelDisplacement_mm_;
FLDamperLength = FrontLeft.DamperLength_mm_;

RRWheelDisplacment = RearRight.WheelDisplacement_mm_;
RRDamperLength = RearRight.DamperLength_mm_;

RLWheelDisplacment = RearLeft.WheelDisplacement_mm_;
RLDamperLength = RearLeft.DamperLength_mm_;


%% Find Quantified Motion Ratio for each corner

MotionRatioFR = GetMotionRatio(FRWheelDisplacment, FRDamperLength);
MotionRatioFL = GetMotionRatio(FLWheelDisplacment, FLDamperLength);
MotionRatioRR = GetMotionRatio(RRWheelDisplacment, RRDamperLength);
MotionRatioRL = GetMotionRatio(RLWheelDisplacment, RLDamperLength);

AMRFR = mean(MotionRatioFR);
AMRFL = mean(MotionRatioFL);
AMRRR = mean(MotionRatioRR);
AMRRL = mean(MotionRatioRL);

% Print motion ratios
disp('--- Motion Ratios ---')
fprintf('Front Left  Motion Ratio: %.3f\n', AMRFL)
fprintf('Front Right Motion Ratio: %.3f\n', AMRFR)
fprintf('Rear Left   Motion Ratio: %.3f\n', AMRRL)
fprintf('Rear Right  Motion Ratio: %.3f\n', AMRRR)


%% Plot Motion Ratio vs Wheel Travel

figure(2)
clf

colors = {'b','r','g','m'};

% ---- Front Left ----
subplot(2,2,1)
plot(FLWheelDisplacment, MotionRatioFL, ...
    'Color', colors{1}, 'LineWidth', 1.5)
xlabel('Wheel Displacement (mm)')
ylabel('Motion Ratio (dDamper / dWheel)')
title('Front Left')
grid on

% ---- Front Right ----
subplot(2,2,2)
plot(FRWheelDisplacment, MotionRatioFR, ...
    'Color', colors{2}, 'LineWidth', 1.5)
xlabel('Wheel Displacement (mm)')
ylabel('Motion Ratio (dDamper / dWheel)')
title('Front Right')
grid on

% ---- Rear Left ----
subplot(2,2,3)
plot(RLWheelDisplacment, MotionRatioRL, ...
    'Color', colors{3}, 'LineWidth', 1.5)
xlabel('Wheel Displacement (mm)')
ylabel('Motion Ratio (dDamper / dWheel)')
title('Rear Left')
grid on

% ---- Rear Right ----
subplot(2,2,4)
plot(RRWheelDisplacment, MotionRatioRR, ...
    'Color', colors{4}, 'LineWidth', 1.5)
xlabel('Wheel Displacement (mm)')
ylabel('Motion Ratio (dDamper / dWheel)')
title('Rear Right')
grid on

sgtitle('Instantaneous Motion Ratio vs Wheel Travel')


%% Functions

function MotionRatio = GetMotionRatio(displacement, DamperLength)

    dDamper = gradient(DamperLength);
    dWheel  = gradient(displacement);

    MotionRatio = dDamper ./ dWheel;

end

