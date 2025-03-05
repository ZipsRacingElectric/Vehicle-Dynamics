%{
This script takes raw downforce and drag data, converts it to matlab data
ready for simulink lookup tables.
%}

% Load the spreadsheet from the specified sheet
filename = 'zr25_aero_data.xlsx';
sheetname = 'matlab';
data = readtable(filename, 'Sheet', sheetname);

% load vehicle data to get CG position
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

% CG dist from axles - these numbers were used to calculate axle loads
% Note: they didn't use ZR25 dimensions to caluclate downforce!! We will
% correct this by calculating a CoP relative to actual ZR25 CG position,
% and then apply those to the actual axle positions in simulink.
front_axle_to_center = 0.762; % center of wheelbase, not CG because we love complication!!!
rear_axle_to_center = 0.762;

% Extract data columns
aero_data_velocity = data{:, 1}; % Velocity (m/s)
front_downforce_no_drs = data{:, 2}; % Front axle downforce (No DRS)
front_downforce_drs = data{:, 3}; % Front axle downforce (With DRS)
rear_downforce_no_drs = data{:, 4}; % Rear axle downforce (No DRS)
rear_downforce_drs = data{:, 5}; % Rear axle downforce (With DRS)
drag_no_drs = data{:, 6}; % Drag (No DRS)
drag_drs = data{:, 7}; % Drag (With DRS)

% Compute total downforce
downforce_no_drs = front_downforce_no_drs + rear_downforce_no_drs;
downforce_drs = front_downforce_drs + rear_downforce_drs;

% Compute center of pressure relative to center of wheelbase
cop_no_drs = (front_downforce_no_drs * front_axle_to_center - rear_downforce_no_drs * rear_axle_to_center) ./ downforce_no_drs;
cop_drs = (front_downforce_drs * front_axle_to_center - rear_downforce_drs * rear_axle_to_center) ./ downforce_drs;

% Calculate offset between center of wheelbase and cg
offset = zr25.a - zr25.wheelbase / 2;

% Calculate center of pressure relative to CG location
cop_no_drs = cop_no_drs + offset;
cop_drs = cop_drs + offset;

% Calculate and offset between center of wheelbase and ZR25 CG position

% Save the results
save('zr25_aero_data.mat', 'aero_data_velocity', 'downforce_no_drs', 'downforce_drs', ...
    'drag_no_drs', 'drag_drs', 'cop_no_drs', 'cop_drs');

% ----- Generate Plots -----

% Plot Downforce and Drag vs. Velocity
figure;
hold on;
plot(aero_data_velocity, downforce_no_drs, '-o', 'DisplayName', 'Downforce (No DRS)');
plot(aero_data_velocity, downforce_drs, '--o', 'DisplayName', 'Downforce (DRS)');
plot(aero_data_velocity, drag_no_drs, '-s', 'DisplayName', 'Drag (No DRS)', 'Color', 'r');
plot(aero_data_velocity, drag_drs, '--s', 'DisplayName', 'Drag (DRS)', 'Color', 'r');
xlabel('Velocity (m/s)');
ylabel('Force (N)');
title('Downforce and Drag vs. Velocity');
legend('Location', 'best');
grid on;
hold off;

% Plot Center of Pressure vs. Velocity
figure;
hold on;
plot(aero_data_velocity, cop_no_drs, '-o', 'DisplayName', 'CoP (No DRS)');
plot(aero_data_velocity, cop_drs, '--o', 'DisplayName', 'CoP (DRS)');
xlabel('Velocity (m/s)');
ylabel('Center of Pressure (m)');
title('Center of Pressure vs. Velocity');
legend('Location', 'best');
grid on;
hold off;
