% Wheel MOI Thing
% takes a velocity vector and estimates the cost per lap of larger tires vs
% heavier tires, based on the rotational kinetic energy

% Written by Logan Haydu April 2025

% filepath = '..\DATA_BUCKET\SampleDataForDataBucket.csv';
Testdata = readmatrix("C:\Users\ATuck\OneDrive - The University of Akron\ZR25\Vehicle Dynamics\Data\2025.06.16 - Mary Gladwin Shakedown\decoded_vd.csv"');  % or whatever your file is called
time = Testdata(:,1);                        

% First column = time
% Average RPM from 4 wheels
AVGRPM = (Testdata(:,2) + Testdata(:,3) + Testdata(:,4) + Testdata(:,5)) / 4;

% Effective radius of the wheel (meters)
r = (13 / 2) * 0.0254;  % 13-inch wheel diameter converted to radius in meters

% Calculate angular velocity (rad/s)
omega = 2 * pi * AVGRPM / 60;

% Calculate linear velocity (m/s)
v_mps = omega * r;
v = v_mps


%% want regen or no want regen
regen_enabled = false;     % Set to false to disable regen
regen_efficiency = 0.6;   % Only matters if regen_enabled = true



%% Time and Velocity:
%time = data.Time;
% v = data.VehicleSpeed / 2.237;

% Column Vector Correction:
time = time(:);
v = v(:);

% Car Parameters:
m = 230;

% Tire Parameters
I1 = .25         % guestimate, change tomorrow ask teams pretty pls tonks
I2 = .46;       % what unit is pls confirm

r_tire1_m = 16 * 25.4/1000/2
r_tire2_m = 20. * 25.4/1000/2

% Translational Kinetic Energy:
KE_t = 0.5 * m *v.^2;

% Rotational Kinetic Energy:
w1 = v / r_tire1_m;
KE_r1_one = 0.5 * I1 * w1.^2;
KE_r1_all = 4 * KE_r1_one;
w2 = v / r_tire2_m;
KE_r2_one = 0.5 * I2 * w2.^2;
KE_r2_all = 4 * KE_r2_one;

%% ENERGY DIFFERENCE CALCULATIONS %%

% Time Step
dt = time(2) - time(1);

% Angular Acceleration
alpha1 = [0; diff(w1)/dt];
alpha2 = [0; diff(w2)/dt];

% Torque and Power
tau1 = I1 * alpha1;
tau2 = I2 * alpha2;
P1 = tau1 .* w1;
P2 = tau2 .* w2;

% Energy accounting
energy_consumed_r1 = 0;
energy_consumed_r2 = 0;
cumulative_energy_r1 = zeros(size(time));
cumulative_energy_r2 = zeros(size(time));

for i = 1:length(time)
    % For Tire 1
    if tau1(i) >= 0  % Acceleration: energy consumed
        energy_consumed_r1 = energy_consumed_r1 + P1(i) * dt;
    else  % Deceleration
        if regen_enabled
            energy_consumed_r1 = energy_consumed_r1 + (1 - regen_efficiency) * abs(P1(i)) * dt;
        else
            energy_consumed_r1 = energy_consumed_r1 + abs(P1(i)) * dt;
        end
    end
    cumulative_energy_r1(i) = energy_consumed_r1;

    % For Tire 2
    if tau2(i) >= 0
        energy_consumed_r2 = energy_consumed_r2 + P2(i) * dt;
    else
        if regen_enabled
            energy_consumed_r2 = energy_consumed_r2 + (1 - regen_efficiency) * abs(P2(i)) * dt;
        else
            energy_consumed_r2 = energy_consumed_r2 + abs(P2(i)) * dt;
        end
    end
    cumulative_energy_r2(i) = energy_consumed_r2;
end

% 4 wheels on car
cumulative_energy_r2 = cumulative_energy_r2 .*4;
cumulative_energy_r1 = cumulative_energy_r1 .*4;


%% RESULTS AND PLOTS %%
disp(['Total Energy Consumed (r = 10" tires): ' num2str(energy_consumed_r1/1000) ' kJ']);
disp(['Total Energy Consumed (r = 13" tires): ' num2str(energy_consumed_r2/1000) ' kJ']);

% Velocity vs Time Plot
figure(1);
plot(time,v, 'b', 'LineWidth', 1);
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
plot(time, .001.*cumulative_energy_r1, 'b-', 'LineWidth', 1); hold on;
plot(time, .001.*cumulative_energy_r2, 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Cumulative Energy Consumed (kJ)');
title('Cumulative Energy Consumption Over 1 lap (with Regen)');
legend('r = 10" Tires', 'r = 13" Tires');
grid on;


totalDist = cumtrapz(time,v);
totalDist = totalDist(end);

scaleyfactorydoo = 22000/totalDist; % total distance on data under endurance distance

differenceeydoo = (cumulative_energy_r2(end) - cumulative_energy_r1(end) ) * scaleyfactorydoo % total kJ diff with regen

differenceeydoo = differenceeydoo/3600 % in kwh