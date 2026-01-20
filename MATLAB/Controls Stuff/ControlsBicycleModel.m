% bicycle model??? 
% Written 01/12/26 Abigail Tucker

clc
clear
close all


%% Bicycle model parameters

m = 313.5;        % kg
Lf = 0.7344;       % m
Lr = 0.7956;       % m
L = Lf+Lr ;

%% Load test data

Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Bicycle Model data test.csv");

% Define the subset of rows to analyze
rows = 19327:21276;

% Extract and smooth the subset
U  = movmean(Data.SPEED(rows), length(rows));
r  = movmean(Data.BOSCH_Z_ANGLE_RATE(rows), length(rows));
ay = movmean(Data.BOSCH_Y_ACCELERATION(rows), length(rows));
time = Data.timestamps(rows);


%% Lateral Load Estimation

[Fyf, Fyr] = GetLateralLoads(m, U, r, Lf, Lr); % newtons

%% Slip angle estimate

[SAF, SAR] = GetSlipAngles(Lf, Lr, U, r);

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
U_mean  = mean(U);        % m/s
ay_mean = mean(ay(~isnan(ay)));   % m/s² (~1.3 g for skidpad)

% --- Understeer gradient ---
K_us = (Lf*Cf_mean - Lr*Cr_mean) / (Cf_mean + Cr_mean);  % meters

% --- Steady-state turn radius from geometry ---
R_turn = U_mean^2 / ay_mean;   % meters

% --- Ideal yaw rate ---
% r_ideal = U_mean / R_turn;         % rad/s
% r_ideal_deg = rad2deg(r_ideal);    % deg/s

% --- Ideal tire steering angle including understeer ---
delta_tire_rad = (L / R_turn) + K_us * (ay_mean / 9.81);  % rad
delta_tire_deg = rad2deg(delta_tire_rad);               % deg

delta_wheel_deg = delta_tire_deg * SteeringRatio;       % deg
SteeringAngle = linspace(0, 80, length(Ca_r))';  % sweep slightly larger than nominal

% Convert steering to road wheel angle (rad)
delta = deg2rad(SteeringAngle) ./ SteeringRatio;

% Denominator (critical-speed sensitive)
den = (Ca_f .* Ca_r .* L^2) - ...
      (m .* U.^2) .* (Lf .* Ca_f - Lr .* Ca_r);

% Ideal yaw rate (rad/s)
r_ideal = delta .* (Ca_f .* Ca_r .* L .* U ./ den);

% Optional: convert to deg/s
r_ideal_deg = rad2deg(r_ideal);


% --- Convert to steering wheel angle for lookup ---
% delta_wheel_deg = delta_tire_deg * SteeringRatio;       % deg

% --- Generate a linear sweep around nominal steering wheel angle ---
% SteeringAngle_deg = linspace(0, delta_wheel_deg*50, 1000);  % sweep slightly larger than nominal
% r_lookup_deg = (SteeringAngle_deg / delta_wheel_deg) .* r_ideal_deg;  % linear interpolation

%% --- Plot lookup ---
figure(2);
plot(SteeringAngle, r_ideal_deg, 'b', 'LineWidth', 1.5); hold on
grid on
xlabel('Steering Wheel Angle [deg]')
ylabel('Ideal Yaw Rate [deg/s]')
title('Steady-State Lookup Table: Ideal Yaw Rate vs Steering Wheel Angle')

% --- Export lookup table ---
lookup_table = table(SteeringAngle_deg', r_lookup_deg', 'VariableNames', {'SteeringWheel_deg','IdealYawRate_deg_s'});

%% Ideal yaw (2) not correct!!
% 
% SteeringRatio = 5.764;            % wheel / tire
% 
% % --- Steering sweep (steering wheel angle) ---
% SteeringWheelAngle = linspace(0, 90, 3001);  
% 
% % Convert to tire angle
% TireSteeringAngle = SteeringWheelAngle./SteeringRatio;
% 
% % find ideal r
% r_Ideal = TireSteeringAngle.*(((Ca_f.*Ca_r).*L.*U)/((Ca_f.*Ca_r.*(L.^2))-((m.*(U.^2)).*((Lf.*Ca_f)-(Lr.*Ca_r)))));
% 
% figure(2)
% plot(SteeringWheelAngle,r_Ideal)
% grid on
% xlabel('Steering Wheel Angle')
% ylabel('Yaw Rate deg/s')

%% Functions

function [LateralLoadFront, LateralLoadRear] = GetLateralLoads(Masskg, Speedkmperhr, YawRatedegpers, LengthFrontm, LengthRearm)
    % This function estimates lateral loads for a bicycle model based on mass,
    % yaw, speed, and vehicle geometry. Fy = (mur)/(1+Lx/Ly)
    
    %   Converts km/h → m/s and deg/s → rad/s
    
    % -------- Unit conversions --------
    U = Speedkmperhr * (1000/3600);         % km/h → m/s
    r = YawRatedegpers * (pi/180);         % deg/s → rad/s
    
    LateralLoadFront = ((Masskg.*U).*r)./(1+(LengthFrontm/LengthRearm));
    LateralLoadRear = ((Masskg.*U).*r)./(1+(LengthRearm/LengthFrontm));
    
    % output in Newtons

end


function [SlipangleFront, SlipAngleRear] = GetSlipAngles(LengthFrontm, LengthRearm, Speedkmperhr, YawRatedegpers)

  % -------- Unit conversions --------
    U = Speedkmperhr * (1000/3600);         % km/h → m/s
    r = YawRatedegpers * (pi/180);         % deg/s → rad/s


    SlipangleFront = (LengthFrontm.*r)./U ;
    SlipAngleRear = (LengthRearm.*r)./U ;

    % Output radians

end