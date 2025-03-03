clear; close all; clc;

% load data from AMK
amk_80C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T80C.mat");
amk_100C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T100C.mat");
amk_120C = load("./DD5 motor characteristic diagram/A2370DD_Matlab/A2370DD_T120C.mat");

fields = fieldnames(amk_80C);
for i = 1:numel(fields)
    field_name = fields{i};
    data = cat(3, amk_80C.(field_name), amk_100C.(field_name), amk_120C.(field_name));
    assignin('base', field_name, data);
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

% Determine max Id from Iq_Peak
max_id = max(Iq_Peak(:));

%% Field Weakening Implementation
v_batt = [300, 400, 500, 600];
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

%% Graph Results (With and Without Field Weakening)

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