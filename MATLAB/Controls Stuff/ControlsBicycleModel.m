% bicycle model??? 
% Written 01/12/26 Abigail Tucker

%% Bicycle model parameters
% not very accurate!!!
m = 250;        % kg
Iz = 142.92;    % kg*m^2
Lf = 0.7;       % m
Lr = 0.83;       % m
L = Lf+Lr ;
%% Cornering stiffnesses
% System of equations
% m(dv+V*r) = Fy1 +Fy2
% 0 = Lf*Fy1-Lr*Fy2



Ca_f = ;         % N/rad  
Ca_r = ;         % N/rad  

%% ========================================================================
%% 1) VELOCITY SWEEP  (small constant steering)
%% ========================================================================
delta = 2 * pi/180;  % 2 degrees steering input (radians)
u_vec = linspace(5,40,25); % m/s sweep (11–90 mph)

yawRate_vel = zeros(size(u_vec));
beta_vel    = zeros(size(u_vec));

for i = 1:length(u_vec)
    u = u_vec(i);

    % Solve steady-state algebraic bicycle equations:
    % beta and r such that dv=0 & dr=0

    % State equations (linear bicycle):
    % 0 = -(Ca_f + Ca_r)/(m*u)*beta + (Ca_r*Lr - Ca_f*Lf)/(m*u)*r + Ca_f/m*delta
    % 0 = (Ca_r*Lr - Ca_f*Lf)/(Iz)*beta - (Ca_f*Lf^2 + Ca_r*Lr^2)/(Iz*u)*r + Ca_f*Lf/Iz*delta

    A = [ -(Ca_f+Ca_r)/(m*u),     (Ca_r*Lr - Ca_f*Lf)/(m*u);
          (Ca_r*Lr - Ca_f*Lf)/Iz, -(Ca_f*Lf^2 + Ca_r*Lr^2)/(Iz*u) ];

    B = [ Ca_f/m*delta;
          Ca_f*Lf/Iz*delta ];

    x = A\B;  % solve [beta; r]
    beta_vel(i)    = x(1);
    yawRate_vel(i) = x(2);
end

%% ========================================================================
%% 2) STEERING SWEEP (fixed velocity)
%% ========================================================================
u = 20;                               % pick one speed (m/s)
delta_vec = linspace(0,5*pi/180,30);  % 0→5 degrees

yawRate_steer = zeros(size(delta_vec));
beta_steer    = zeros(size(delta_vec));

for i = 1:length(delta_vec)
    delta_i = delta_vec(i);

    A = [ -(Ca_f+Ca_r)/(m*u),     (Ca_r*Lr - Ca_f*Lf)/(m*u);
          (Ca_r*Lr - Ca_f*Lf)/Iz, -(Ca_f*Lf^2 + Ca_r*Lr^2)/(Iz*u) ];

    B = [ Ca_f/m*delta_i;
          Ca_f*Lf/Iz*delta_i ];

    x = A\B;  % solve [beta; r]
    beta_steer(i)    = x(1);
    yawRate_steer(i) = x(2);
end

%% ========================================================================
%% PLOTTING
%% ========================================================================

figure;
plot(u_vec, yawRate_vel, 'LineWidth',2)
xlabel('Velocity (m/s)')
ylabel('Yaw Rate (rad/s)')
title('Velocity Sweep — Yaw Rate vs Speed')
grid on

figure;
plot(delta_vec*180/pi, yawRate_steer, 'LineWidth',2)
xlabel('Steering Angle (deg)')
ylabel('Yaw Rate (rad/s)')
title('Steering Sweep — Yaw Rate vs Steering Input (u = 20 m/s)')
grid on

disp('Done! You now have yaw rate vs speed AND yaw rate vs steering.')
