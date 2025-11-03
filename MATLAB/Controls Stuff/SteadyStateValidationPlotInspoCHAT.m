%% 4WD EV Steady-State Validation Sweep
% Sweeps steering angles and speeds, computes lateral acceleration and yaw rate
% Vehicle: your parameter object

% --- Define sweep ranges ---
delta_deg = -10:1:10;      % steering angles in degrees
Vx_sweep  = [5, 10, 15];   % longitudinal speeds [m/s]

% Preallocate results
ay_map = zeros(length(Vx_sweep), length(delta_deg));
r_map  = zeros(length(Vx_sweep), length(delta_deg));

% Loop over speeds and steering angles
for i = 1:length(Vx_sweep)
    Vx = Vx_sweep(i);
    for j = 1:length(delta_deg)
        delta = delta_deg(j) * pi/180;  % convert to radians

        % Call steady-state solver
        [Vy, r, ay] = steadyStateVehicleLevel(Vehicle, delta, Vx);

        % Store results
        ay_map(i,j) = ay;
        r_map(i,j)  = r;
    end
end

%% --- Plot Results ---
figure;
subplot(2,1,1);
plot(delta_deg, ay_map', 'LineWidth', 1.5);
xlabel('Steering Angle δ [deg]');
ylabel('Lateral Acceleration a_y [m/s^2]');
title('Steady-State Lateral Acceleration vs Steering Angle');
legend(arrayfun(@(v) ['V_x = ', num2str(v), ' m/s'], Vx_sweep, 'UniformOutput', false));
grid on;

subplot(2,1,2);
plot(delta_deg, r_map', 'LineWidth', 1.5);
xlabel('Steering Angle δ [deg]');
ylabel('Yaw Rate r [rad/s]');
title('Steady-State Yaw Rate vs Steering Angle');
legend(arrayfun(@(v) ['V_x = ', num2str(v), ' m/s'], Vx_sweep, 'UniformOutput', false));
grid on;
