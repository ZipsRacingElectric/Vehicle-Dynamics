clear, clc

run('/Users/khoi/Documents/GitHub/Vehicle-Dynamics/MATLAB/Kinematics+compliance_tests/MotionRatioDataProcessor.m');
%Given Motion Ratio, Spring Rate

ZR25 = vehicle('/Users/khoi/Documents/GitHub')

springRate = [0.4, 0.5, 0.1, 10, 40, 25, 70, 150, 500, 1200]

frontTrackWidth = 50
readTrackWidth = 55

wheelRateFR = springRate .* MotionRatioFR.^2
wheelRateFL = springRate .* MotionRatioFL.^2
wheelRateRR = springRate .* MotionRatioRR.^2
wheelRateRL = springRate .* MotionRatioRL.^2

