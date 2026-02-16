% bicycle model 🚲🚲🚲
% Written 01/12/26 Abigail Tucker

clc
clear
close all


%% Bicycle model parameters

% load in vehicle object 🏎️
addpath vehicle_data ;
githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr25_data.xlsx');
ZR25 = vehicle(parameterSpreadsheet);

m = ZR25.mass_total;        % kg
a = ZR25.a;       % m
b = ZR25.b;       % m
L = ZR25.wheelbase ;

%% Load test data

Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Bicycle Model data test.csv");
Data = Data(Data.timestamps >= 130 & Data.timestamps <= 150, :); % define time range of interest

%% data + filtering

time = Data.timestamps; % s

%Time / Sampling 
dt = mean(diff(time));
Fs = 1/dt;                 % Sampling frequency (Hz)

% Raw Signals 
U_raw  = Data.SPEED;                    % km/h
r_raw  = Data.BOSCH_Z_ANGLE_RATE;       % deg/s
ay_raw = Data.BOSCH_Y_ACCELERATION;     % g

% Unit Conversions 
U_raw_con  = U_raw * (1000/3600);   % km/h -> m/s
r_raw_con  = r_raw * (pi/180);      % deg/s -> rad/s
ay_raw_con = ay_raw * 9.81;         % g -> m/s^2

%  Low-Pass Butterworth Filter
fc = 6;                         % Cutoff frequency (Hz)
[bfilt,afilt] = butter(2, fc/(Fs/2));    % 2nd-order Butterworth

% Apply Zero-Phase Filtering 
U  = filtfilt(bfilt,afilt,U_raw_con);
r  = filtfilt(bfilt,afilt,r_raw_con);
ay = filtfilt(bfilt,afilt,ay_raw_con);


%% Lateral Load Estimation

[Fyf, Fyr] = GetLateralLoads(m, ay, a, b); % newtons

%% Slip angle estimate
delta_est = (L .* r) ./ U;   % rad

[SAF, SAR] = GetSlipAngles(a, b, U, r, delta_est); % rad

SAF_deg = rad2deg(SAF);   %deg negative included for convention, need to look into more
SAR_deg = rad2deg(SAR);   %deg

%% body slip estimate

beta = zeros(size(U));
tau = 5;                 % decay time constant (seconds)
alpha = dt/tau;          % decay factor

for i = 2:length(U)

    if U(i) > 5
        betadot_i = (ay(i)/U(i)) - r(i);
    else
        betadot_i = 0;
    end

    % Leaky integration
    beta(i) = beta(i-1) + betadot_i*dt - alpha*beta(i-1);

end

beta_deg = rad2deg(beta);


%% Cornering stiffnesses
% System of equations
% m(dv+V*r) = Fy1 +Fy2
% 0 = Lf*Fy1-Lr*Fy2

Ca_f = Fyf./SAF;         % N/rad  +
Ca_r = Fyr./SAR;         % N/rad  

Ca_f_deg = (Ca_f)* (pi/180); %N/deg
Ca_r_deg = (Ca_r)* (pi/180); %N/deg

%% Plot cornering stiffness and speed

figure(1);

% --- Cornering stiffness ---
subplot(2,1,1)
plot(time, Ca_f_deg, 'b', 'LineWidth', 1.2)
hold on
plot(time, Ca_r_deg, 'r', 'LineWidth', 1.2)

% Mean lines
yline_f = mean(Ca_f_deg(~isnan(Ca_f_deg)));  % ignore NaNs
yline_r = mean(Ca_r_deg(~isnan(Ca_r_deg)));

yline(yline_f, '--k', ['Front Tire 🏎️ Mean = ' num2str(round(yline_f)) ' N/deg'], 'LabelHorizontalAlignment','left','FontSize',10)
yline(yline_r, '--k', ['Rear Tire  🏎️ Mean = ' num2str(round(yline_r)) ' N/deg'], 'LabelHorizontalAlignment','left','FontSize',10)

grid on
xlabel('Time [s]')
ylabel('Cornering Stiffness C_\alpha [N/deg]')
legend('Front', 'Rear')
title('Front and Rear Cornering Stiffness')

% --- Vehicle speed ---
subplot(2,1,2)
plot(time, U, 'g', 'LineWidth',1.2)
grid on
xlabel('Time [s]')
ylabel('Speed [m/s]')
title('Vehicle Speed')
hold off

%% Steady-State Lookup Table: Steering Wheel Angle → Ideal Yaw Rate (1)

% --- Vehicle parameters ---
Cf_mean = mean(Ca_f(~isnan(Ca_f)));  
Cr_mean = mean(Ca_r(~isnan(Ca_r)));  

SteeringRatio = 5.764;   % column-to-tire ratio

% --- Steady-state mean speed and lateral acceleration ---
%U_mean  = mean(U);  %m/s
ay_mean = mean(ay(~isnan(ay)));   % m/s^2


% Understeer gradient (rad/(m/s^2))
K_us = (m * (b*Cr_mean - a*Cf_mean)) / (Cf_mean * Cr_mean * L);



% Ideal Steering Angle Sweep
aygen = linspace(0, 2*9.81, 200);   % 0 to 2g sweep (m/s^2)

U_ref = mean(U(~isnan(U)));          % representative speed (scalar)

delta = (L/(U_ref^2) + K_us) .* aygen;   % rad
delta_deg = rad2deg(delta);               % deg

% Numerator denominator 
Num = delta*(Cf_mean*Cr_mean*L.*U_ref); % 
Den = (Cf_mean*Cr_mean*(L^2))-(m*(U_ref^2)*(a*Cf_mean-b*Cr_mean));

r_ideal = Num/Den; % rad/s
r_ideal_deg = r_ideal*(180/pi);

SteeringWheelAngle = delta_deg*SteeringRatio;


%% Plot lookup 🔎
figure(2);
plot(SteeringWheelAngle, r_ideal_deg, 'b', 'LineWidth', 1.5); hold on
grid on
xlabel('Steering Wheel Angle [deg]')
ylabel('Ideal Yaw Rate [deg/s]')
title('Steady-State: Ideal Yaw Rate vs Steering Wheel Angle')

% --- Export lookup table ---
lookup_table = table(SteeringWheelAngle', r_ideal_deg', 'VariableNames', {'SteeringWheel_deg','IdealYawRate_deg_s'});


%% Plot tire 🛞 
figure(3);
scatter(SAF_deg, Fyf, 'blue', 'LineWidth', 1.5);
hold on
scatter(SAR_deg, Fyr, 'red','LineWidth',1.5);
xlim([-10 10]);
ylim([-1500 1500]);
hold on
grid on
xlabel('slip angle')
ylabel('lateral force')
title('lateral force vs slip angle');


figure(4);
plot(time, ay_raw, 'blue', 'LineWidth', 1.5);
hold on
grid on
xlabel('time')
ylabel('lateral accel')
title('lateral accel vs time');
%% Cubic fit tire model 🧊

valid_f = ~isnan(SAF_deg) & ~isnan(Fyf);
valid_r = ~isnan(SAR_deg) & ~isnan(Fyr);

% Fit cubic using only valid data
Af = [SAF_deg(valid_f), SAF_deg(valid_f).^3];
Ar = [SAR_deg(valid_r), SAR_deg(valid_r).^3];

coef_f = Af \ Fyf(valid_f);
coef_r = Ar \ Fyr(valid_r);

C1_f = coef_f(1);
C3_f = -coef_f(2);

C1_r = coef_r(1);
C3_r = -coef_r(2);

% ----- Limit model to measured data range -----
max_alpha_f = max(abs(SAF_deg(valid_f)));
max_alpha_r = max(abs(SAR_deg(valid_r)));

alpha_plot_f = linspace(-max_alpha_f, max_alpha_f, 200);
alpha_plot_r = linspace(-max_alpha_r, max_alpha_r, 200);

Fy_front_model = CubicTireModel(alpha_plot_f, C1_f, C3_f);
Fy_rear_model  = CubicTireModel(alpha_plot_r, C1_r, C3_r);


% ----- Plot -----
figure(5);
scatter(SAF_deg(valid_f), Fyf(valid_f), 10, 'b'); hold on
scatter(SAR_deg(valid_r), Fyr(valid_r), 10, 'r');

plot(alpha_plot_f, Fy_front_model, 'k','LineWidth',2);
plot(alpha_plot_r, Fy_rear_model, 'k--','LineWidth',2);


xlabel('Slip Angle [deg]')
ylabel('Lateral Force [N]')
title('Measured Data vs Cubic Tire Model')
grid on


%% sanity plots
figure(6)
plot(time, beta_deg)
grid on
xlabel('time')
ylabel('deg')

%% Functions

function [LateralLoadFront, LateralLoadRear] = GetLateralLoads(Masskg, ay, a, b)
    % This function estimates lateral loads for a bicycle model.
    % input mass - kg, ay - m/s^2, lf, lr - meters

    TotalLatForce = Masskg*ay; % Netwtons

    LateralLoadFront = (b*TotalLatForce)/(a+b);
    LateralLoadRear = (a*TotalLatForce)/(a+b);
    
    % output in Newtons

end

function [SlipAngleFront, SlipAngleRear] = GetSlipAngles(a, b, speed, r, delta)

  
    % Minimum speed threshold (m/s)
    minSpeed = 5;   % adjust if needed (~11 mph)

    % Preallocate with NaN (so bad regions don't explode)
    SlipAngleFront = NaN(size(speed));
    SlipAngleRear  = NaN(size(speed));

    % Valid data mask
    valid = speed > minSpeed;

    % Slip angles (rad)
    SlipAngleFront(valid) = (a .* r(valid)) ./ speed(valid) - delta(valid);
    SlipAngleRear(valid)  = -(b .* r(valid)) ./ speed(valid);

end

function Fy = CubicTireModel(alpha, C1, C3)
    Fy = C1 .* alpha - C3 .* alpha.^3;
end
