%% Motion Ratio 

% this code will plot vertical wheel travel against corresponding damper
% length providing the motion ratio for each corner 

% Written B=by Abigail Tucker 01/09/26

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


%% Plot Corners

figure(1)

% Define colors and markers for each curve
colors = {'b','r','g','m'};

subplot(2, 2, 1)
plot(FLWheelDisplacment, FLDamperLength,'Color', colors{1}, 'LineWidth',1.5)
xlabel('Damper Length (mm)'); ylabel('Wheel Displacement, (mm)')
grid on

subplot(2, 2, 2)
plot(FRWheelDisplacment, FRDamperLength, 'Color', colors{2}, 'LineWidth',1.5)
xlabel('Damper Length (mm)'); ylabel('Wheel Displacement, (mm)')
grid on

subplot(2, 2, 3)
plot( RLWheelDisplacment, RLDamperLength, 'Color', colors{3}, 'LineWidth',1.5)
xlabel('Damper Length (mm)'); ylabel('Wheel Displacement, (mm)')
grid on

subplot(2, 2, 4)
plot( RRWheelDisplacment, RRDamperLength, 'Color', colors{4}, 'LineWidth',1.5)
xlabel('Damper Length (mm)'); ylabel('Wheel Displacement, (mm)')
grid on

% Add a single title above all subplots
sgtitle('Motion Ratio By Corner')

%% Functions

function  MotionRatio = GetMotionRatio(displacement, DamperLength)

            MotionRatio = displacement./DamperLength;
end 

