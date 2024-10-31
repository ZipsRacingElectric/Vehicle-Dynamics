clc; clear all;

wheelbase = 1.53;
track_width = linspace(1, 1.3, 100);
contact_patch_deflection = linspace(0, 0.08, 100)';

%% Condition: wheels are straight
% For Fy there is no error. For TV moments where we are applying Fx:
tire_long_force = 1; % magnitude of force has no effect on error

expected_yaw_moment = tire_long_force * (track_width/2);

% columns: increasing track width. rows: increasing contact patch deflection
actual_yaw_moment = tire_long_force * ((track_width/2) - contact_patch_deflection);


yaw_moment_error = (actual_yaw_moment - expected_yaw_moment) ./ expected_yaw_moment .* 100;


figure;
surf(track_width, contact_patch_deflection, yaw_moment_error)
title('Yaw Moment Error vs Long. Force vs Contact Patch Deflection')
xlabel('Track Width m')
ylabel('Contact Patch Deflection m')
zlabel('Yaw Moment Error %')
%}