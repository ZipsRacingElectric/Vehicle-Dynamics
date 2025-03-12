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
Note: all of these tables are functions of current and speed
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

% create dimension vectors
motor_temps = [80, 100, 120];
motor_currents = linspace(0, 105, 21);
motor_speeds = linspace(0, 20000, 201);

% Create interpolation functions
torque_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Shaft_Torque);
voltage_phase_rms_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Voltage_Phase_RMS);
vq_rms_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Vq_RMS);
vd_rms_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Vd_RMS);
loss_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Total_Loss);
Iq_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Iq_Peak);
Id_interp = griddedInterpolant({motor_speeds, motor_currents, motor_temps}, Id_Peak);

% Determine max Id from Iq_Peak
max_id = max(Iq_Peak(:));

%% Field Weakening Implementation
v_batt = [100, 200, 300, 400, 500, 600];
motor_temperature = 80;
torque = zeros(length(v_batt), length(motor_speeds));
power = zeros(length(v_batt), length(motor_speeds));

for i = 1:length(v_batt)
    phase_voltage_limit = v_batt(i) / sqrt(2);
    
    for u = 1:length(motor_speeds)
        speed = motor_speeds(u);
        
        % Iterate over possible Id values
        max_torque = 0;
        best_iq = 0;
        best_id = 0;
        
        id_step = max_id / 20;
        id_test = 0;
        while id_test <= max_id % Sweep field weakening current
            id_test = id_test + id_step;
            
            % Calculate total phase voltage with Id
            iq_test = sqrt(max(motor_currents.^2 - id_test^2, 0));
            vq_test = vq_rms_interp(speed, max(iq_test(:)), motor_temperature);
            vd_test = vd_rms_interp(speed, max(id_test(:)), motor_temperature);
            v_phase_test = sqrt(vd_test.^2 + vq_test.^2);
            
            % Ensure voltage limit is not exceeded
            if v_phase_test <= phase_voltage_limit
                torque_test = torque_interp(speed, max(iq_test(:)), motor_temperature);
                if torque_test > max_torque
                    max_torque = torque_test;
                    best_iq = iq_test;
                    best_id = id_test;
                end
            end
        end
        
        % Store results
        torque(i, u) = max_torque;
        power(i, u) = max_torque * (speed * 2 * pi / 60);
    end
end

%% Compute Efficiency Map for Each DC Bus Voltage
for u = 1:length(v_batt)
    figure;
    efficiency_map = NaN(length(motor_speeds), length(motor_speeds)); % Initialize efficiency map
    
    for i = 1:length(motor_speeds)
        for j = 1:length(torque(u, :))
            % Get speed and torque values
            speed = motor_speeds(i);
            torque_value = min(torque(u, j), 40);
            
            % Ensure within feasible operating range
            if torque_value > torque(u, i) || torque_value < 0
                continue;
            end
            
            % Compute Mechanical Power
            P_mech = torque_value * (speed * 2 * pi / 60);
            
            % Compute Electrical Power Input
            Iq = Iq_interp(speed, torque_value, motor_temps(1));
            Id = Id_interp(speed, torque_value, motor_temps(1));
            Vq = vq_rms_interp(speed, Iq, motor_temps(1));
            Vd = vd_rms_interp(speed, Id, motor_temps(1));
            P_loss = loss_interp(speed, Iq, motor_temps(1));
            
            % Compute Efficiency
            efficiency_map(i, j) = P_mech / (P_mech + P_loss);
        end
    end
    
    % Plot Efficiency Contour
    [X, Y] = meshgrid(motor_speeds, linspace(0, max(torque(:)), length(motor_speeds)));
    contour(X, Y, efficiency_map, 20, 'LineWidth', 1.5);
    colorbar;
    xlabel("Motor Speed (RPM)");
    ylabel("Motor Torque (Nm)");
    title(sprintf("Motor Efficiency Contour - DC Bus Voltage %dV", v_batt(u)));
    caxis([0 100]);
    grid on;
end

%% Graph Results
% Plot Torque vs Speed without Field Weakening
figure;
hold on;
no_fw_torque = zeros(length(v_batt), length(motor_speeds));
for i = 1:length(v_batt)
    phase_voltage = v_batt(i) / sqrt(2);
    for u = 1:length(motor_speeds)
        speed = motor_speeds(u);
        voltage_slice = squeeze(voltage_phase_rms_interp((zeros(1, length(motor_currents)) + speed), motor_currents, (zeros(1, length(motor_currents)) + motor_temperature)));
        tolerance = 0.1;
        valid_idx = abs(voltage_slice - phase_voltage) < tolerance;
        if any(valid_idx)
            currents = max(motor_currents(valid_idx));
        elseif all(voltage_slice < phase_voltage)
            currents = motor_currents(end);
        elseif all(voltage_slice > phase_voltage)
            currents = 0;
        else
            currents = interp1(voltage_slice, motor_currents, phase_voltage);
        end
        no_fw_torque(i, u) = torque_interp(speed, currents, motor_temperature);
    end
    plot(motor_speeds, no_fw_torque(i, :), '--', 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV (No FW)', v_batt(i)));
end
xlabel("Motor Speed (RPM)"); ylabel("Motor Torque (Nm)");
title("Torque vs Speed without Field Weakening");
grid on; legend; hold off;


% Plot Torque vs Speed without Field Weakening
figure;
hold on;
no_fw_power = zeros(length(v_batt), length(motor_speeds));
for i = 1:length(v_batt)
    for u = 1:length(motor_speeds)
        no_fw_power(i, u) = no_fw_torque(i, u) * (motor_speeds(u) * 2 * pi / 60);
    end
    plot(motor_speeds, no_fw_power(i, :), '--', 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV (No FW)', v_batt(i)));
end
xlabel("Motor Speed (RPM)"); ylabel("Motor Power (W)");
title("Power vs Speed without Field Weakening");
grid on; legend; hold off;


% Plot Torque vs Speed without Field Weakening
figure;
hold on;
for i = 1:length(v_batt)
    plot(motor_speeds, torque(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV', v_batt(i)));
end
xlabel("Motor Speed (RPM)"); ylabel("Motor Torque (Nm)");
title("Torque vs Speed with Field Weakening");
grid on; legend; hold off;

% Plot Torque vs Speed without Field Weakening
figure;
hold on;
for i = 1:length(v_batt)
    plot(motor_speeds, power(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Bus Voltage = %dV', v_batt(i)));
end
xlabel("Motor Speed (RPM)"); ylabel("Motor Power (W)");
title("Power vs Speed with Field Weakening");
grid on; legend; hold off;