%% Written by Carl Fahrner
%%Bumpsteer code, calculating bumpsteer gain using rideheight travel and toe angle 
%% Estamates ride height sensitivity 
clear, clc
%% 
Excel=readtable("C:\Users\Fahrn\OneDrive\Desktop\Ecar team\Bumpsteer excel.xlsx") ;
%% 
WheelLift=Excel.ZMovement_mm_ ;
%% 
D=Excel.DistanceOnWall ;
%% 
Dial=Excel.DialIndicator ;
%% 
Distancetowall=Excel.DistanceToWall_m_ ;
%% 
y=Distancetowall .* 1000 ;
%% 
BumpsteerAngle=atand((D - Dial) ./ (y)) ;
%% 
figure(1)
subplot (3,1,1);
plot(WheelLift,BumpsteerAngle);
xlabel("ZMovement_mm_");
ylabel('Bumpsteer Results (degrees)');
title('Bumpsteer Results vs Distance to Wall');
grid on;
%% 
BumpsteerGain = gradient(BumpsteerAngle, WheelLift);

subplot(3,1,2);
plot(WheelLift, BumpsteerGain);
xlabel("ZMovement_mm_");
ylabel('Bumpsteer Gain (degrees/mm)');
title('Bumpsteer Gain vs ZMovement');
grid on;

yl=ylim;

operatingRange = patch([30 80 80 30],...
[yl(1) yl(1) yl(2) yl(2)],...
[0.85 0.85 0.85],...
'FaceAlpha', 0.5, ...
'EdgeColor', 'b');

RideHight = 60;

xline(RideHight, '--k', 'RideHeight', 'Linewidth', 1.5)

operatingwindow = [30 80];
opIdx = WheelLift >= operatingwindow(1) & WheelLift <= operatingwindow(2);

meanGain = mean(BumpsteerGain(opIdx));
stdGain = std(BumpsteerGain(opIdx));

subplot(3,1,3);
plot(WheelLift(opIdx), BumpsteerGain(opIdx))
yline(meanGain, ':b', 'meanGain', 'LineWidth', 1.0)
xlabel("ZMovement_mm_");
ylabel('Bumpsteer Gain (degrees/mm)');
ylim([-0.01,0.01])
title('Bumpsteer Gain vs ZMovement');
grid on;
%% 
fprintf('meanGain: %.4f\n', meanGain)
fprintf('Ride height tuning sensitivity: %.4f\n', stdGain)