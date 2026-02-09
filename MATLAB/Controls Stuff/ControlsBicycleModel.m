% bicycle model??? 
% Written 01/12/26 Abigail Tucker

clc
clear
close all


%% Bicycle model parameters

addpath vehicle_data ;
githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr26_data.xlsx');
ZR26 = vehicle(parameterSpreadsheet);

m = 313.5;        % kg
Lf = 0.7344;       % m
Lr = 0.7956;       % m
L = Lf+Lr ;

%% Load test data

Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Bicycle Model data test.csv");
Data = Data(Data.timestamps >= 175 & Data.timestamps <= 258, :);

% Define the subset of rows to analyze
% Extract and smooth the subset
U = (Data.SPEED);
r = (Data.BOSCH_Z_ANGLE_RATE);
ay = (Data.BOSCH_Y_ACCELERATION);

time = Data.timestamps; % s


%% Lateral Load Estimation

[Fyf, Fyr] = GetLateralLoads(m, ay, Lf, Lr); % newtons

%% Slip angle estimate

[SAF, SAR] = GetSlipAngles(Lf, Lr, U, r); % rad

SAF_deg = -rad2deg(SAF);   %deg
SAR_deg = -rad2deg(SAR);   %deg

%% Cornering stiffnesses
% System of equations
% m(dv+V*r) = Fy1 +Fy2
% 0 = Lf*Fy1-Lr*Fy2

Ca_f = Fyf./SAF;         % N/rad  
Ca_r = Fyr./SAR;         % N/rad  

%% Plot cornering stiffness and speed

figure(1);

% --- Cornering stiffness ---
subplot(2,1,1)
plot(time, Ca_f, 'b', 'LineWidth', 1.2)
hold on
plot(time, Ca_r, 'r', 'LineWidth', 1.2)

% Mean lines
yline_f = mean(Ca_f(~isnan(Ca_f)));  % ignore NaNs
yline_r = mean(Ca_r(~isnan(Ca_r)));

yline(yline_f, '--k', ['Front Tire 🏎️ Mean = ' num2str(round(yline_f)) ' N/rad'], 'LabelHorizontalAlignment','left','FontSize',10)
yline(yline_r, '--k', ['Rear Tire  🏎️ Mean = ' num2str(round(yline_r)) ' N/rad'], 'LabelHorizontalAlignment','left','FontSize',10)

grid on
xlabel('Time [s]')
ylabel('Cornering Stiffness C_\alpha [N/rad]')
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
U_mean  = mean(U);        % km/hr
ay_mean = mean(ay(~isnan(ay)));   % g

U_meanm = U_mean*(1000/3600); % km/h -> m/s
ay_meanm = ay_mean*(9.81); % g -> m/s


% Understeer gradient (rad/(m/s^2))
K_us = (m * (Lr*Cr_mean - Lf*Cf_mean)) / (Cf_mean * Cr_mean * L);



% Ideal Steering Angle Sweep
aygen = linspace(0, 2*9.81, 200);   % 0 to 2g sweep (m/s)

delta = (L/U_meanm^2 + K_us) .* aygen; % rad 
delta_deg = delta*(180/pi); % convert to deg


% Numerator denominator 
Num = delta*(Cf_mean*Cr_mean*L*U_meanm); % 
Den = (Cf_mean*Cr_mean*(L^2))-(m*(U_meanm^2)*(Lf*Cf_mean-Lr*Cr_mean));

r_ideal = Num/Den; % rad/s
r_ideal_deg = r_ideal*(180/pi);

SteeringWheelAngle = delta_deg*SteeringRatio;


%% --- Plot lookup ---
figure(2);
plot(SteeringWheelAngle, r_ideal_deg, 'b', 'LineWidth', 1.5); hold on
grid on
xlabel('Steering Wheel Angle [deg]')
ylabel('Ideal Yaw Rate [deg/s]')
title('Steady-State: Ideal Yaw Rate vs Steering Wheel Angle')

% --- Export lookup table ---
lookup_table = table(SteeringWheelAngle', r_ideal_deg', 'VariableNames', {'SteeringWheel_deg','IdealYawRate_deg_s'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Plot lookup ---
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
plot(time, ay, 'blue', 'LineWidth', 1.5);
hold on
grid on
xlabel('time')
ylabel('lateral accel')
title('lateral accel vs time');
%% Cubic fit tire model 

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
figure;
scatter(SAF_deg(valid_f), Fyf(valid_f), 10, 'b'); hold on
scatter(SAR_deg(valid_r), Fyr(valid_r), 10, 'r');

plot(alpha_plot_f, Fy_front_model, 'k','LineWidth',2);
plot(alpha_plot_r, Fy_rear_model, 'k--','LineWidth',2);


xlabel('Slip Angle [deg]')
ylabel('Lateral Force [N]')
title('Measured Data vs Cubic Tire Model')
grid on




%% Functions

function [LateralLoadFront, LateralLoadRear] = GetLateralLoads(Masskg, ayg, Lfm, Lrm)
    % This function estimates lateral loads for a bicycle model.
 
       
    Ay = ayg*9.81;

    TotalLatForce = Masskg*Ay; % Netwtons

    LateralLoadFront = (Lrm*TotalLatForce)/(Lfm+Lrm);
    LateralLoadRear = (Lfm*TotalLatForce)/(Lfm+Lrm);
    
    % output in Newtons

end

function [SlipAngleFront, SlipAngleRear] = GetSlipAngles(LengthFrontm, LengthRearm, Speedkmperhr, YawRatedegpers)

    % Convert units
    U = Speedkmperhr * (1000/3600);   % km/h → m/s
    r = YawRatedegpers * (pi/180);    % deg/s → rad/s

    % Minimum speed threshold (m/s)
    minSpeed = 5;   % adjust if needed (~11 mph)

    % Preallocate with NaN (so bad regions don't explode)
    SlipAngleFront = NaN(size(U));
    SlipAngleRear  = NaN(size(U));

    % Valid data mask
    valid = U > minSpeed;

    % Slip angles (rad)
    SlipAngleFront(valid) = (LengthFrontm .* r(valid)) ./ U(valid);
    SlipAngleRear(valid)  = (LengthRearm  .* r(valid)) ./ U(valid);

end

function Fy = CubicTireModel(alpha, C1, C3)
    Fy = C1 .* alpha - C3 .* alpha.^3;
end
