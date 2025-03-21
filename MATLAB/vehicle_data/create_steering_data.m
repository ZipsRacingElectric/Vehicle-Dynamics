%{
This script takes spreadsheet data from OptimumLap and creates workspace
arrays which can be loaded into simulink lookup tables
%}

% Load the spreadsheet from the specified sheet
filename = 'zr25_steering_data.xlsx';
sheetname = 'matlab';
data = readtable(filename, 'Sheet', sheetname);

% Extract data columns
steering_data_sw = data{:, 1}; % Velocity (m/s)
steering_data_fl = data{:, 2}; % Front axle downforce (No DRS)
steering_data_fr = data{:, 3}; % Front axle downforce (With DRS)

% Save the results
save('zr25_steering_data.mat', 'steering_data_sw', 'steering_data_fl', 'steering_data_fr');

% ----- Generate Plots -----
figure;
hold on;
plot(steering_data_sw, steering_data_fl, 'DisplayName', 'FL Steer Angle (deg)');
plot(steering_data_sw, steering_data_fr, 'DisplayName', 'FR Steer Angle (deg)');
xlabel('Steering Wheel Angle (deg)');
ylabel('Tire Angles (deg)');
title('Steering System Curves ZR25');
legend('Location', 'best');
grid on;
hold off;
