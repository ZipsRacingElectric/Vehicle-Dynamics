% bicycle model
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

% need to fix units and convert after loading in

Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Bicycle Model data test.csv");

% Define the subset of rows to analyze
rows = 19327:21276;

% Extract and smooth the subset
U = (Data.SPEED(rows));
r = (Data.BOSCH_Z_ANGLE_RATE(rows));
ay = (Data.BOSCH_Y_ACCELERATION(rows));


fs = 100;      % sampling frequency

fc = 5;        % cutoff frequency (Hz)

[b,a] = butter(2, fc/(fs/2));   % 2nd order Butterworth

ay = filtfilt(b,a, ay);
r  = filtfilt(b,a, r);
U  = filtfilt(b,a, U);



time = Data.timestamps(rows); % s

plot(time, Data.BOSCH_Y_ACCELERATION)

%% Lateral Load Estimation

[Fyf, Fyr] = GetLateralLoads(m, ay, Lf, Lr); % newtons




%% Slip angle estimate

[SAF, SAR] = GetSlipAngles(Lf, Lr, U, r); % rad

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

r_mean = mean(r(~isnan(r)));   % g


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


%% Body Slip Estimation 
r_mean = movmean(r, length(time));  % moving average over full vector


BetaDot_rad = (ay_mean.*9.81)./(U_mean.*(1000/3600))-(r_mean*(pi/180));
BetaDot_deg = BetaDot_rad * (180/pi);
BetaDot_deg = BetaDot_deg - mean(BetaDot_deg);


dt = 1/100;   % 100 Hz
beta_deg = cumtrapz(BetaDot_deg) * dt;

r_rad = r_mean * pi/180;

% Body slip at front/rear (degrees)
beta_f_deg = beta_deg + atan2d(Lf * r_rad, U_mean*(1000/3600));
beta_r_deg = beta_deg - atan2d(Lr * r_rad, U_mean*(1000/3600));


%% --- Plot lookup ---
figure(2);
plot(SteeringWheelAngle, r_ideal_deg, 'b', 'LineWidth', 1.5); hold on
grid on
xlabel('Steering Wheel Angle [deg]')
ylabel('Ideal Yaw Rate [deg/s]')
title('Steady-State: Ideal Yaw Rate vs Steering Wheel Angle')

% --- Export lookup table ---
lookup_table = table(SteeringWheelAngle', r_ideal_deg', 'VariableNames', {'SteeringWheel_deg','IdealYawRate_deg_s'});


figure(3)
subplot(4,1,1)
plot(time, beta_f_deg)
title ('beta')

subplot(4,1,2)
plot(time, ay)
title ('ay')


subplot(4,1,3)
plot(time, r)
title ('r')

subplot(4,1,4)
plot(time, U)
title ('u')
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
 
    U = Speedkmperhr * (1000/3600);            % km/h → m/s
    r = YawRatedegpers * (pi/180);             % deg/s → rad/s

    % Slip angles (rad) 
    SlipAngleFront = (LengthFrontm .* r) ./ U;  % rad
    SlipAngleRear  = (LengthRearm  .* r) ./ U;  % rad
end
