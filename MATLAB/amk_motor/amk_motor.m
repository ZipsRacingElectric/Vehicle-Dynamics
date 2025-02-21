%{

%% Overview
This script takes the supplied AMK test data, and develops speed - torque
curves and calculates maximum motor torque for a given battery voltage and
motor temperature.

%% Details:
- V_peak = V_rms * sqrt(2). DC bus voltage is what limits V_peak.
Therefore, for a given DC bus voltage V_dc, V_rms_line = V_dc / sqrt(2)

%% TODO:
- graph motor losses and calculate electrical power
- graph torque/power curves with field weakening enabled
- graph efficiency map of motor
- save data for simulink lookup table. Re-calculate current as
function of fixed phase voltage and save as matrix so we don't need to do
any math in simulink

%% Data Description:
- Speed: motor shaft speed [RPM]
- Shaft_Torque: motor shaft torque [Nm]
- Stator_Current_Phase_Peak: Amplituide of phase current [A]
- Stator_Current_Phase_RMS: RMS value of phase current [A]
- Stator_Current_Line_Peak: Amplituide of line current [A]
- Stator_Current_Line_RMS: RMS value of line current [A]
- Voltage_Phase_Peak: Peak voltage of motor phase [V]
- Voltage_Phase_RMS: RMS voltage of motor phase [V]
- Voltage_Line_Peak: Peak of line voltage [V]
- Voltage_Line_RMS: RMS of line voltage [V]
- Id_Peak: Amplituide of field weakening current [A]
- Id_RMS: RMS of field weakening current [A]
- Iq_Peak: Amplituide of torque generating current [A]
- Iq_RMS: RMS of torque generating current [A]
- Vd_Peak: Amplituide of field weakening voltage [V]
- Vd_RMS: RMS of field weakening voltage [V]
- Vq_Peak: Amplituide of torque generating voltage [V]
- Vq_RMS: RMS of torque generating voltage [V]
- Frequency: Electrical frequency [Hz]
- Total Loss: Total motor losses [W]
- Stator_Copper_Loss: Resistive losses [W]
- Iron_Loss: Loss from  excitation producing eddy currents in iron [W]
- Magnet_Loss: Magnetic losses [W]
- Mechanical_Loss: frictional, vibration losses [W]
- Power_Factor: [unitless]
- Electromagnetic_Torque: internal torque before subration of losses [Nm]

%}

clear; close all; clc;

% load data from AMK
amk_80C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T80C.mat");
amk_100C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T100C.mat");
amk_120C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T120C.mat");

% Create 3D matrix of data with dimensions ordered as: motor_speeds, motor_currents,
% motor_temps
fields = fieldnames(amk_80C);
for i = 1:numel(fields)
    field_name = fields{i}; % Get field name
    data = cat(3, amk_80C.(field_name), amk_100C.(field_name), amk_120C.(field_name));  % Concatenate field data from all elements
    assignin('base', field_name, data); % Assign to base workspace
end

clear amk_100C amk_120C amk_80C i data fields field_name

% create dimension vectors for breakpoints
motor_temps = [80, 100, 120];
motor_currents = linspace(0, 105, 21);
motor_speeds = linspace(0, 20000, 201);

% Create interpolation functions
torque_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Shaft_Torque);
voltage_phase_rms_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Voltage_Phase_RMS);

%% Save data for simulink lookup table

%% Calculate torque-speed curve for different DC bus voltages, no field weakening, fixed motor temp
v_batt = [300, 400, 500, 600];
motor_temperature = 80;

torque = zeros(length(v_batt), length(motor_speeds));
power = zeros(length(v_batt), length(motor_speeds));

for i = 1:length(v_batt)

    % Calculate motor voltage from DC bus voltage
    phase_voltage = v_batt(i) / sqrt(2);
    
    % For each motor speed, solve for the max current at the fixed motor
    % voltage
    currents = zeros(1, length(motor_speeds));

    for u = 1:length(motor_speeds)
        speed = motor_speeds(u);

        % Extract the 1D row of voltage vs. current at a fixed temperature
        voltage_slice = squeeze(voltage_phase_rms_interp((zeros(1, length(motor_currents)) + speed), motor_currents, (zeros(1, length(motor_currents)) + motor_temperature))); 
        
        % Identify where the voltage matches the fixed voltage within tolerance
        tolerance = 0.1; % Small tolerance for floating-point comparisons
        valid_idx = abs(voltage_slice - phase_voltage) < tolerance;

        if any(valid_idx)
            % this handles edge case where the voltage flat-lines for
            % increasing current and interpolation gets funky

            % Select the highest current where voltage is within tolerance
            currents(u) = max(motor_currents(valid_idx));

        elseif all(voltage_slice < phase_voltage)
            % in this condition, the maximum torque generating current is
            % availible at this speed
            currents(u) = motor_currents(end);

        elseif all(voltage_slice > phase_voltage)
            % in this condition, the motor can not generate any current
            % because the phase voltage is too low for this speed
            currents(u) = 0;

        else
            % Interpolate the current corresponding to fixed voltage
            currents(u) = interp1(voltage_slice, motor_currents, phase_voltage);
        end
    end

    % find torque generated for maximum currents
    torque(i, :) = torque_interp(motor_speeds, currents, (zeros(1, length(motor_speeds)) + motor_temperature));

    % calulate mechanical power
    power(i, :) = torque(i, :) .* (motor_speeds .* 2.*pi() ./ 60);

end

%% Graph Results
figure;
hold on;

for i = 1:length(v_batt)
    plot(motor_speeds, torque(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV', v_batt(i)));
end

% Labels & Title
xlabel("Motor Speed (RPM)", 'FontSize', 12, 'FontWeight', 'bold');
ylabel("Motor Torque (Nm)", 'FontSize', 12, 'FontWeight', 'bold');
title("Torque vs Speed Curve - No Field Weakening, T = 80 deg C", 'FontSize', 14, 'FontWeight', 'bold');
xticks = get(gca, 'XTick');
yticks = get(gca, 'YTick');
xticklabels(arrayfun(@(x) sprintf('%d', x), xticks, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%d', y), yticks, 'UniformOutput', false));
xlim([0 22000]);   
ylim([0 30]);   
grid on;
grid minor;
legend('Location', 'northeast', 'FontSize', 10);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 12;
ax.GridAlpha = 0.2;  % Adjust transparency (lower = lighter)
ax.MinorGridAlpha = 0.1; % Even lighter minor grid
hold off;

figure;
hold on;

for i = 1:length(v_batt)
    plot(motor_speeds, power(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV', v_batt(i)));
end

% Labels & Title
xlabel("Motor Speed (RPM)", 'FontSize', 12, 'FontWeight', 'bold');
ylabel("Motor Power (W)", 'FontSize', 12, 'FontWeight', 'bold');
title("Power vs Speed Curve - No Field Weakening, T = 80 deg C", 'FontSize', 14, 'FontWeight', 'bold');
xticks = get(gca, 'XTick');
yticks = get(gca, 'YTick');
xticklabels(arrayfun(@(x) sprintf('%d', x), xticks, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%d', y), yticks, 'UniformOutput', false));
xlim([0 22000]);   
ylim([0 44000]);   
grid on;
grid minor;
legend('Location', 'northeast', 'FontSize', 10);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 12;
ax.GridAlpha = 0.2;  % Adjust transparency (lower = lighter)
ax.MinorGridAlpha = 0.1; % Even lighter minor grid
hold off;
