% Wheel MOI Thing
% takes a velocity vector and estimates the cost per lap of larger tires vs
% heavier tires, based on the rotational kinetic energy

% Written by Logan Haydu April 2024

filepath = '..\DATA_BUCKET\SampleDataForDataBucket.csv'
columnTitles = string({'Time', 'VehicleSpeed'})

data = Data_Brita(filepath,columnTitles)



% Time and Velocity:
time = data.Time


% Car Parameters:
m = 200;

% Tire Parameters:
m_tire1 = 10;
r_tire1_in = 10;
r_tire1_m = r_tire1_in / 39.37;
I1 = 0.5 * m_tire1 * r_tire1_m^2;
m_tire2 = 13;
r_tire2_in = 13;
r_tire2_m = r_tire2_in / 39.37;
I2 = 0.5 * m_tire2 * r_tire2_m^2;

% Translational Kinetic Energy:
KE_t = 0.5 * m *data.VehicleSpeed.^2;

% Rotational Kinetic Energy:
w1 = data.VehicleSpeed / r_tire1_m;
KE_r1_one = 0.5 * I1 * w1.^2;
KE_r1_all = 4 * KE_r1_one;
w2 = data.VehicleSpeed / r_tire2_m;
KE_r2_one = 0.5 * I2 * w2.^2;
KE_r2_all = 4 * KE_r2_one;

%% ENERGY DIFFERENCE CALCULATIONS %%

% Time Step
dt = time(2) - time(1);

% Angular Acceleration
alpha1 = [0 diff(w1)/dt];
alpha2 = [0 diff(w2)/dt];

% Torque and Power
tau1 = I1 * alpha1;
tau2 = I2 * alpha2;
P1 = tau1 .* w1;
P2 = tau2 .* w2;

% Energy accounting
energy_consumed_r1 = 0;
energy_consumed_r2 = 0;
regen_efficiency = 0.6;  % 60% recovery
regen_enabled = true;
cumulative_energy_r1 = zeros(size(time));
cumulative_energy_r2 = zeros(size(time));

for i = 1:length(time)
    % Tire 1 (Acceleration or Deceleration)
    if tau1(i) >= 0  % If accelerating
        energy_consumed_r1 = energy_consumed_r1 + P1(i) * dt;
    else  % Decelerating (Regenerating)
        energy_consumed_r1 = energy_consumed_r1 + (1 - regen_efficiency * regen_enabled) * abs(P1(i)) * dt;  % Save energy during regen
    end
    cumulative_energy_r1(i) = energy_consumed_r1;
    
    % Tire 2 (Acceleration or Deceleration)
    if tau2(i) >= 0  % If accelerating
        energy_consumed_r2 = energy_consumed_r2 + P2(i) * dt;
    else  % Decelerating (Regenerating)
        energy_consumed_r2 = energy_consumed_r2 + (1 - regen_efficiency * regen_enabled) * abs(P2(i)) * dt;  % Save energy during regen
    end
    cumulative_energy_r2(i) = energy_consumed_r2;
end

%% RESULTS AND PLOTS %%
disp(['Total Energy Consumed (r = 10" tires): ' num2str(energy_consumed_r1) ' J']);
disp(['Total Energy Consumed (r = 13" tires): ' num2str(energy_consumed_r2) ' J']);

% Velocity vs Time Plot
figure(1);
plot(time,data.VehicleSpeed, 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time');
grid on;

% Translational Kinetic Energy Plot
figure(2);
subplot(2,1,1);
plot(time, KE_t, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy (J)');
title('Translational Kinetic Energy (Car)');
legend('Translational KE');
grid on;

% Rotational Kinetic Energy Comparisons Plot
subplot(2,1,2);
plot(time, KE_r1_all, 'g--', 'LineWidth', 1); hold on;
plot(time, KE_r1_one, 'k--', 'LineWidth', 1);
plot(time, KE_r2_all, 'm-', 'LineWidth', 1);
plot(time, KE_r2_one, 'c-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy (J)');
title('Rotational Kinetic Energy (Based on Tire Size)');
legend('4 Tires, r = 10"', '1 Tire, r = 10"', ...
       '4 Tires, r = 13"', '1 Tire, r = 13"');
grid on;

% Cumulative Energy Plot
figure(3);
plot(time, cumulative_energy_r1, 'b-', 'LineWidth', 1); hold on;
plot(time, cumulative_energy_r2, 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Cumulative Energy Consumed (J)');
title('Cumulative Energy Consumption Over Time (with Regen)');
legend('r = 10" Tires', 'r = 14" Tires');
grid on;