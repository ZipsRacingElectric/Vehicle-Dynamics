%{

%% Overview
This script creates and plots math models representing compliance test data

%% Details:
All compliances effect the overall vehicle damping coefficient which should
be accuretly modeled because that is what torque vectoring active controllers
 are trying to improve upon. See the compliance manifestor ppts in the 
VD research folder on the OneDrive. The following two are the biggest
 contributors, because changes in slip angle (anyting toe related) result
 in the largest changes in the tire forces, and because these compliances
 individually act on only one side of the car, resulting in a balance
 change.

 Mz Front Steer Compliance (aka ENF): important because deflection due to bad 
steering rack/column design. Generally non-linear. The data is from 
Vishnu Sanjay from 4ze racing, it should be of the same magnituide until 
we can measure our own cars.

Fy Rear Steer Compliance (aka EYR): biggest contribution to rear toe
compliance. Generally linear. Singificant factor in transient response.
Typical values are 0.1 deg/1000 N for FSAE (again, should be measured).

Rear Roll Steer Compliance: not modeled, could have big influence on rear
suspension.

Fx Rear Steer Compliance: not modeled - deflection of toe due to accel/braking of rear
wheels. Also important factor, however not modeled because eliminating bump
steer will have a bigger effect on toe reduction from Fx forces.

%% TODO:
- quantify these compliances with actual test data

%}

%% ENF Compliance
% this is a non-linear model for ENF compliance. It is common to use a
% simple linear model and specify compliance as [deg per 100 N-m]

% 1st column: Mz (N-m), 2nd column: toe deflection (deg)
enf_data = csvread("mz_steer_compliance_test_data.csv");

% fit a polynomial to the data
enf_model = polyfit(enf_data(:,1), enf_data(:,2), 5);

% view raw data and model
mz_grid = linspace(-100, 100, 100);
enf_deflection = polyval(enf_model, mz_grid);

figure;
hold on;
plot(enf_data(:,1), enf_data(:,2));
plot(mz_grid, enf_deflection);
title("Front toe compliance due to Mz moment");
xlabel("Mz (N-m)");
ylabel("Toe Deflection, per wheel (deg)");
hold off;

%% EYR Compliance
% Note: depending on how the toe link is designed on the car, Fy force will
% either contribute to (+) or (-) toe change resulting in an over or
% under steering effect. Typically (+) means toe out and (-) is toe-in
eyr_gain = 0.1; % [deg per 1000 N]

fy_grid = linspace(-2500, 2500, 100);
fyr_deflection = fy_grid * eyr_gain/1000;

figure;
plot(fy_grid, fyr_deflection);
title("Rear toe compliance due to Fy force");
xlabel("Fy (N)");
ylabel("Toe Deflection, per wheel (deg)");

