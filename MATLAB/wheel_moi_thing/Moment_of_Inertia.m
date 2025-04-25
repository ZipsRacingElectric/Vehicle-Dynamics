% Wheel MOI Thing
% takes a velocity vector and estimates the cost per lap of larger tires vs
% heavier tires, based on the rotational kinetic energy blah blah

% written by Logan Haydu April 2024


% Time and Velocity:
time = 0:100;
v = (3/4)*cos((time/4) + pi) + 0.75;

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
KE_t = 0.5 * m * v.^2;

% Rotational Kinetic Energy:
w1 = v / r_tire1_m;
KE_r1_one = 0.5 * I1 * w1.^2;
KE_r1_all = 4 * KE_r1_one;
w2 = v / r_tire2_m;
KE_r2_one = 0.5 * I2 * w2.^2;
KE_r2_all = 4 * KE_r2_one;

% Velocity vs Time Plot
figure(1);
plot(time, v, 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time');
grid on;

% Translational Kinetic Energy Plot
figure(2);
subplot(2,1,1);
plot(time, KE_t, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Energy (Joules)');
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
ylabel('Energy (Joules)');
title('Rotational Kinetic Energy (Based on Tire Size)');
legend('4 Tires, r = 10"', '1 Tire, r = 10"', ...
       '4 Tires, r = 14"', '1 Tire, r = 14"');
grid on;