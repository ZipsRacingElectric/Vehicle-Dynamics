% bicycle model 🚲🚲🚲
% Written 01/12/26 Abigail Tucker

% all of this could really be tucked into one "getbicyclemmodel" function
% that can be called to any time spread determined manually in assmadif.
% Only meant to be accurate for steady state in generally linear region of
% tires. Should not need to be run live on car, dump lookup table into
% simulink 

% FOR FUTURE USERS: you need to update vehicle data spreadsheet, steering
% ratio (hidden in function below).Youll need to clear my 4 circles and
% adapt to whatever your chosen data set consists of (this is where one
% funciton would be way easier). 

% Could be expanded to double track,, tucked into one function, better
% filtering and interpolation methods, lots of other stuff.Cant decide if
% interpolated real steering ang or a simulated sweep makes more sense. 

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
Iz = ZR25.yaw_polar_inertia; %kg * m^2,

%% Load test data

Data = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Data Excels\Goodyear data slalom + circles.csv");
Data = Data(Data.timestamps >= 0 & Data.timestamps <= 500, :); % define time range of interest
%time = Data.timestamps; % s

%% data + filtering

% define circles
smallCirc = Data(Data.timestamps >= 25 & Data.timestamps <= 31, :); % data seems mehhhh
MediumCirc = Data(Data.timestamps >= 61.5 & Data.timestamps <= 75, :);
LargeCirc = Data(Data.timestamps >= 86 & Data.timestamps <= 104, :);
XLargeCirc = Data(Data.timestamps >= 120 & Data.timestamps <= 142, :);



% Filter data by section
[U, r, ay, time] = GetGoodData(Data); % full set
[U_Small, r_Small, ay_Small, time_Small] = GetGoodData(smallCirc);
[U_Medium, r_Medium, ay_Medium, time_Medium] = GetGoodData(MediumCirc);
[U_Large, r_Large, ay_Large, time_Large] = GetGoodData(LargeCirc);
[U_XLarge, r_XLarge, ay_XLarge, time_XLarge] = GetGoodData(XLargeCirc);


%% Lateral Load Estimation

[Fyf, Fyr] = GetLateralLoads(m, ay, a, b); % newtons
[Fyf_Small, Fyr_Small] = GetLateralLoads(m, ay_Small, a, b); % newtons
[Fyf_Medium, Fyr_Medium] = GetLateralLoads(m, ay_Medium, a, b); % newtons
[Fyf_Large, Fyr_Large] = GetLateralLoads(m, ay_Large, a, b); % newtons
[Fyf_XLarge, Fyr_XLarge] = GetLateralLoads(m, ay_XLarge, a, b); % newtons

%% Slip angle estimate


[SAF, SAR, SAF_deg, SAR_deg] = GetSlipAngles(a, b, U, r, L); % rad
[SAF_Small, SAR_Small, SAF_Deg_Small, SAR_Deg_Small] = GetSlipAngles(a, b, U_Small, r_Small, L); % rad There is a probelm with this dataset
[SAF_Medium, SAR_Medium, SAF_Deg_Medium, SAR_Deg_Medium] = GetSlipAngles(a, b, U_Medium, r_Medium, L); % rad
[SAF_Large, SAR_Large, SAF_Deg_Large, SAR_Deg_Large] = GetSlipAngles(a, b, U_Large, r_Large, L); % rad
[SAF_XLarge, SAR_XLarge, SAF_Deg_XLarge, SAR_Deg_XLarge] = GetSlipAngles(a, b, U_XLarge, r_XLarge, L); % rad



%% Cornering stiffnesses
% System of equations
% m(dv+V*r) = Fy1 +Fy2
% 0 = Lf*Fy1-Lr*Fy2

[Ca_f, Ca_r, Ca_f_deg, Ca_r_deg] = GetCorneringStiffness(Fyf, Fyr, SAF, SAR);
[Ca_f_Small, Ca_r_Small, Ca_f_deg_Small, Ca_r_deg_Small] = GetCorneringStiffness(Fyf_Small, Fyr_Small, SAF_Small, SAR_Small);
[Ca_f_Medium, Ca_r_Medium, Ca_f_deg_Medium, Ca_r_deg_Medium] = GetCorneringStiffness(Fyf_Medium, Fyr_Medium, SAF_Medium, SAR_Medium);
[Ca_f_Large, Ca_r_Large, Ca_f_deg_Large, Ca_r_deg_Large] = GetCorneringStiffness(Fyf_Large, Fyr_Large, SAF_Large, SAR_Large);
[Ca_f_XLarge, Ca_r_XLarge, Ca_f_deg_XLarge, Ca_r_deg_XLarge] = GetCorneringStiffness(Fyf_XLarge, Fyr_XLarge, SAF_XLarge, SAR_XLarge);


%% body slip estimate
%im not gonna deal with this right now
% R = r*(180/pi)./U; % deg/s m/s -> deg/m
% 
% Wr = ZR25.rear_mass;
% alphar = Wr.*(U.^2)./(Ca_r_deg.*R);
% 
% 
% beta = 57.3*(b./R)-alphar;

%% PLot Slip angle and Lateral Force

figure()
tiledlayout(2,2)

ax1 = nexttile;
GetSlipAnglePlots(ax1, SAF_Deg_Small, SAR_Deg_Small, Fyf_Small, Fyr_Small, 'Small Circle')

ax2 = nexttile;
GetSlipAnglePlots(ax2, SAF_Deg_Medium, SAR_Deg_Medium, Fyf_Medium, Fyr_Medium, 'Medium Circle')

ax3 = nexttile;
GetSlipAnglePlots(ax3, SAF_Deg_Large, SAR_Deg_Large, Fyf_Large, Fyr_Large, 'Large Circle')

ax4 = nexttile;
GetSlipAnglePlots(ax4, SAF_Deg_XLarge, SAR_Deg_XLarge, Fyf_XLarge, Fyr_XLarge, 'XL Circle')


sgtitle('Slip Angle vs Lateral Force Across Circle Tests')

%% Plot cornering stiffness and speed

figure()

ax1 = subplot(2,2,1);
GetCorneringStiffnessPlots(ax1, time_Small, Ca_f_deg_Small, Ca_r_deg_Small, 'Small Circle')

ax2 = subplot(2,2,2);
GetCorneringStiffnessPlots(ax2, time_Medium, Ca_f_deg_Medium, Ca_r_deg_Medium, 'Medium Circle')

ax3 = subplot(2,2,3);
GetCorneringStiffnessPlots(ax3, time_Large, Ca_f_deg_Large, Ca_r_deg_Large, 'Large Circle')

ax4 = subplot(2,2,4);
GetCorneringStiffnessPlots(ax4, time_XLarge, Ca_f_deg_XLarge, Ca_r_deg_XLarge, 'XL Circle')

% % --- Vehicle speed ---
% subplot(2,1,2)
% plot(time, U, 'g', 'LineWidth',1.2)
% grid on
% xlabel('Time [s]')
% ylabel('Speed [m/s]')
% title('Vehicle Speed')
% hold off

%% Steady-State Lookup Table: Steering Wheel Angle → Ideal Yaw Rate (1)


[SteeringWheelAngle, r_ideal, r_ideal_deg] = GetIdealYaw(U, L, SAF_deg, Fyf, SAR_deg, Fyr, ay, m, a ,b);
[SteeringWheelAngle_Small, r_ideal_Small, r_ideal_deg_Small] = GetIdealYaw(U_Small, L, SAF_Deg_Small, Fyf_Small, SAR_Deg_Small,Fyr_Small,ay_Small, m, a, b);
[SteeringWheelAngle_Medium, r_ideal_Medium, r_ideal_deg_Medium] = GetIdealYaw(U_Medium, L, SAF_Deg_Medium, Fyf_Medium, SAR_Deg_Medium, Fyr_Medium,ay_Medium, m, a, b);
[SteeringWheelAngle_Large, r_ideal_Large, r_ideal_deg_Large] = GetIdealYaw(U_Large, L, SAF_Deg_Large, Fyf_Large, SAR_Deg_Large, Fyr_Large, ay_Large, m, a, b);
[SteeringWheelAngle_XLarge, r_ideal_XLarge, r_ideal_deg_XLarge] = GetIdealYaw(U_XLarge, L, SAF_Deg_XLarge, Fyf_XLarge, SAR_Deg_XLarge, Fyr_XLarge, ay_XLarge, m, a, b);


% mirror data
% Small circle
SteeringWheelAngle_Small = [SteeringWheelAngle_Small; -SteeringWheelAngle_Small];
r_ideal_Small           = [r_ideal_Small; -r_ideal_Small];
r_ideal_deg_Small       = [r_ideal_deg_Small; -r_ideal_deg_Small];
U_Small                 = [U_Small; U_Small];

% Medium circle
SteeringWheelAngle_Medium = [SteeringWheelAngle_Medium; -SteeringWheelAngle_Medium];
r_ideal_Medium           = [r_ideal_Medium; -r_ideal_Medium];
r_ideal_deg_Medium       = [r_ideal_deg_Medium; -r_ideal_deg_Medium];
U_Medium                 = [U_Medium; U_Medium];

% Large circle
SteeringWheelAngle_Large = [SteeringWheelAngle_Large; -SteeringWheelAngle_Large];
r_ideal_Large           = [r_ideal_Large; -r_ideal_Large];
r_ideal_deg_Large       = [r_ideal_deg_Large; -r_ideal_deg_Large];
U_Large                 = [U_Large; U_Large];

% XLarge circle
SteeringWheelAngle_XLarge = [SteeringWheelAngle_XLarge; -SteeringWheelAngle_XLarge];
r_ideal_XLarge           = [r_ideal_XLarge; -r_ideal_XLarge];
r_ideal_deg_XLarge       = [r_ideal_deg_XLarge; -r_ideal_deg_XLarge];
U_XLarge                 = [U_XLarge; U_XLarge];

%% Understeer Gradient 
% 
% p = polyfit(U(valid), K(valid), 2); % quadratic fit
% K_fit = polyval(p, U_bp);
% 
% figure
% scatter(U(valid), K(valid), 10, 'filled')
% grid on
% xlabel('Speed [m/s]')
% ylabel('Understeer Gradient K')
% title('Understeer Gradient vs Speed')
%% Alternate yaw for sanity 

% 
% r_ideal_check = (U_ref ./ (L + K_us*U_ref^2)) .* delta;
% r_ideal_deg_check = rad2deg(r_ideal_check);

fit = goodnessOfFit(Ca_r,Ca_f,'NRMSE') ;
% 
% Yaw_table_Check = zeros(length(U_bp), length(delta_bp));
% 
% for i = 1:length(U_bp)
%     for j = 1:length(delta_bp)
% 
%         U = U_bp(i);
%         delta = deg2rad(delta_bp(j) / SteeringRatio);
% 
%         K_use = K_mean; % or interp from fit
% 
%         r = (U / (L + K_use * U^2)) * delta;
% 
%         Yaw_table_Check(i,j) = rad2deg(r);
%     end
% end

%% Plot Ideal Yaw 🔎
figure
tiledlayout(2,2)

ax1 = nexttile;
GetYawPlots(ax1, SteeringWheelAngle_Medium, r_ideal_deg_Medium, 'Small Circle')

ax2 = nexttile;
GetYawPlots(ax2, SteeringWheelAngle_Medium, r_ideal_deg_Medium, 'Medium Circle')

ax3 = nexttile;
GetYawPlots(ax3, SteeringWheelAngle_Large, r_ideal_deg_Large, 'Large Circle')

ax4 = nexttile;
GetYawPlots(ax4, SteeringWheelAngle_XLarge, r_ideal_deg_XLarge, 'XL Circle')

sgtitle('Ideal Yaw Rate vs Steering Wheel Angle Across Circle Tests')

%% Create 2-d table

% Combine all datasets
U_all     = [U_Small; U_Medium; U_Large; U_XLarge];
delta_all = [SteeringWheelAngle_Small;
             SteeringWheelAngle_Medium;
             SteeringWheelAngle_Large;
             SteeringWheelAngle_XLarge];
yaw_all = [r_ideal_deg_Small;
           r_ideal_deg_Medium;
           r_ideal_deg_Large;
           r_ideal_deg_XLarge];
% Mirror
delta_all = [delta_all; -delta_all];
yaw_all   = [yaw_all; -yaw_all];
U_all     = [U_all; U_all];


% PLot raw
figure
hold on
grid on

scatter3(SteeringWheelAngle_Small,  U_Small,  r_ideal_deg_Small,  20,'r','filled')
scatter3(SteeringWheelAngle_Medium, U_Medium, r_ideal_deg_Medium, 20,'g','filled')
scatter3(SteeringWheelAngle_Large,  U_Large,  r_ideal_deg_Large,  20,'b','filled')
scatter3(SteeringWheelAngle_XLarge, U_XLarge, r_ideal_deg_XLarge, 20,'k','filled')

xlabel('Steering Wheel Angle [deg]')
ylabel('Speed [m/s]')
zlabel('Ideal Yaw Rate [deg/s]')

legend('Small Circle','Medium Circle','Large Circle','XLarge Circle')

title('Raw Ideal Yaw Data by Circle')

view(45,25)

% Filter bad
valid = ~isnan(U_all) & ~isnan(delta_all) & ~isnan(yaw_all);

U_all = U_all(valid);
delta_all = delta_all(valid);
yaw_all = yaw_all(valid);



U_bp = linspace(min(U_all), max(U_all), 15);      % speed breakpoints
delta_bp = linspace(-100, 100, 40);
[DELTA, U] = meshgrid(delta_bp, U_bp);

F = scatteredInterpolant(delta_all, U_all, yaw_all,'natural','linear');
Yaw_table = F(DELTA, U);


%3-dPlot
% figure
% surf(delta_bp, U_bp, Yaw_table)
% 
% xlabel('Steering Wheel Angle [deg]')
% ylabel('Speed [m/s]')
% zlabel('Ideal Yaw Rate [deg/s]')
% 
% title('Yaw Rate Lookup Table')
% % colorbar
% % shading interp
% 
% hold on
%% Emoji Plot
figure
hold on
grid on

xlabel('Steering Wheel Angle [deg]')
ylabel('Speed [m/s]')
zlabel('Ideal Yaw Rate [deg/s]')
title('Yaw Rate Lookup Table')

view(45,25)

% Set axis limits so text becomes visible
axis([min(delta_bp) max(delta_bp) min(U_bp) max(U_bp) min(Yaw_table(:)) max(Yaw_table(:))])

for i = 1:length(U_bp)
    for j = 1:length(delta_bp)

        z = Yaw_table(i,j);

        if ~isnan(z)
            text(delta_bp(j), U_bp(i), z, '🐸', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize',30)
        end

    end
end
%% Simulink
% Save for Simulink
% Steering_bp = delta_bp(:);
% Speed_bp = U_bp(:);
% Yaw_table = Yaw_table;
% save('YawLookupTable.mat','Steering_bp','Speed_bp','Yaw_table')


%% Sanity Check

% figure
% scatter(r_Small, r_ideal_deg_Small, 15, 'filled')
% hold on
% 
% x = linspace(min(r_measured_deg), max(r_measured_deg));
% plot(x,x,'k--','LineWidth',2)
% 
% grid on
% xlabel('Measured Yaw Rate [deg/s]')
% ylabel('Ideal Yaw Rate [deg/s]')
% title('Measured vs Ideal Yaw Rate')
%% Plot tire 🛞 
figure();
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


figure();
plot(time, ay, 'blue', 'LineWidth', 1.5);
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
figure();
scatter(SAF_deg(valid_f), Fyf(valid_f), 10, 'b'); hold on
scatter(SAR_deg(valid_r), Fyr(valid_r), 10, 'r');

plot(alpha_plot_f, Fy_front_model, 'k','LineWidth',2);
plot(alpha_plot_r, Fy_rear_model, 'k--','LineWidth',2);


xlabel('Slip Angle [deg]')
ylabel('Lateral Force [N]')
title('Measured Data vs Cubic Tire Model')
grid on


%% sanity plots
% figure(6)
% plot(SteeringWheelAngle, r_ideal_deg_check)
% grid on
% xlabel('steering ang')
% ylabel('r deg')

disp([min(r_ideal_deg_Medium) max(r_ideal_deg_Medium)])
disp(sum(isnan(r_ideal_deg_Small)))
disp([min(SteeringWheelAngle_Large) max(SteeringWheelAngle_Large)])

%% Functions

function [LateralLoadFront, LateralLoadRear] = GetLateralLoads(Masskg, ay, a, b)
    % This function estimates lateral loads for a bicycle model.
    % input mass - kg, ay - m/s^2, lf, lr - meters

    TotalLatForce = Masskg*ay; % Netwtons

    LateralLoadFront = (b*TotalLatForce)/(a+b);
    LateralLoadRear = (a*TotalLatForce)/(a+b);
    
    % output in Newtons

end

function [SlipAngleFront, SlipAngleRear, SAF_deg, SAR_deg] = GetSlipAngles(a, b, speed, r, L)

    delta = (L .* r) ./ speed;   % rad
  
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

    SAF_deg = rad2deg(SlipAngleFront);   %deg 
    SAR_deg = rad2deg(SlipAngleRear);   %deg


end

function Fy = CubicTireModel(alpha, C1, C3)
    Fy = C1 .* alpha - C3 .* alpha.^3;
end

function [U, r, ay, time] = GetGoodData(Data)

% this functions takes raw data and cleans it and converts it to standard
% units
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

end

function [Ca_f, Ca_r, Ca_f_deg, Ca_r_deg] = GetCorneringStiffness(Fyf, Fyr, SAF, SAR)

Ca_f = Fyf./SAF;         % N/rad  +
Ca_r = Fyr./SAR;         % N/rad  

Ca_f_deg = (Ca_f)* (pi/180); %N/deg
Ca_r_deg = (Ca_r)* (pi/180); %N/deg

minSlip = deg2rad(0.5);

valid_f = abs(SAF) > minSlip;
valid_r = abs(SAR) > minSlip;

Ca_f = NaN(size(SAF));
Ca_r = NaN(size(SAR));

Ca_f(valid_f) = Fyf(valid_f)./SAF(valid_f);
Ca_r(valid_r) = Fyr(valid_r)./SAR(valid_r);

end

function GetCorneringStiffnessPlots(ax, time, Ca_f_deg, Ca_r_deg, titleStr)

axes(ax)  % tell MATLAB where to plot

plot(time, Ca_f_deg, 'b', 'LineWidth', 1.2)
hold on
plot(time, Ca_r_deg, 'r', 'LineWidth', 1.2)

% Mean lines
yline_f = mean(Ca_f_deg(~isnan(Ca_f_deg)));
yline_r = mean(Ca_r_deg(~isnan(Ca_r_deg)));

yline(yline_f, '--k', ['Front Tire 🏎️ Mean = ' num2str(round(yline_f)) ' N/deg'], 'LabelHorizontalAlignment','right','FontSize',14) 
yline(yline_r, '--k', ['Rear Tire 🏎️ Mean = ' num2str(round(yline_r)) ' N/deg'], 'LabelHorizontalAlignment','right','FontSize',14)

grid on
xlabel('Time [s]')
ylabel('C_\alpha [N/deg]')
title(titleStr)
legend('Front','Rear','Location','best')

end

function GetSlipAnglePlots(ax, SAF_deg, SAR_deg, Fyf, Fyr, titleStr)

axes(ax)
hold on

% Remove NaNs
validF = ~isnan(SAF_deg) & ~isnan(Fyf);
validR = ~isnan(SAR_deg) & ~isnan(Fyr);

SAF_deg = SAF_deg(validF);
Fyf     = Fyf(validF);

SAR_deg = SAR_deg(validR);
Fyr     = Fyr(validR);

% Scatter data
scatter(SAF_deg, Fyf, 12, 'b', 'filled')
scatter(SAR_deg, Fyr, 12, 'r', 'filled')

% Linear fit (cornering stiffness)
pF = polyfit(SAF_deg, Fyf, 1);
pR = polyfit(SAR_deg, Fyr, 1);

xF = linspace(min(SAF_deg), max(SAF_deg),100);
xR = linspace(min(SAR_deg), max(SAR_deg),100);

plot(xF, polyval(pF,xF),'k-','LineWidth',2)
plot(xR, polyval(pR,xR),'k-','LineWidth',2)

grid on
xlabel('Slip Angle [deg]')
ylabel('Lateral Force [N]')
title(titleStr)

legend('Front Tire','Rear Tire','Front Fit','Rear Fit','Location','best')

end

function [SteeringWheelAngle, r_ideal, r_ideal_deg] = GetIdealYaw(U, L, SAF_deg, Fyf, SAR_deg, Fyr, ay, m, a, b)

Cf = polyfit(SAF_deg, Fyf,1);
Cr = polyfit(SAR_deg, Fyr,1);

Cf = Cf(1);
Cr = Cr(1);
SteeringRatio = 5.764;

% Steering from bicycle kinematics
delta = (L .* ay) ./ (U.^2);     % rad
SteeringWheelAngle = rad2deg(delta) * SteeringRatio;

% Ideal yaw
Num = delta .* (Cf * Cr * L .* U);
Den = (Cf * Cr * (L^2)) - (m * (U.^2) * (a*Cf - b*Cr));

r_ideal = Num ./ Den;
r_ideal_deg = rad2deg(r_ideal);

end

function GetYawPlots(ax, SteeringWheelAngle, r_ideal_deg, titleStr)

axes(ax)
hold on

%Remove NaNs
valid = ~isnan(SteeringWheelAngle) & ~isnan(r_ideal_deg) & isfinite(r_ideal_deg);
SteeringWheelAngle = SteeringWheelAngle(valid);
r_ideal_deg = r_ideal_deg(valid);

% Sort (important for plotting)
[SteeringWheelAngle, idx] = sort(SteeringWheelAngle);
r_ideal_deg = r_ideal_deg(idx);

% Plot
plot(SteeringWheelAngle, r_ideal_deg, 'LineWidth', 2, 'Color', 'r')

grid on
xlabel('Steering Wheel Angle [deg]')
ylabel('Ideal Yaw Rate [deg/s]')
title(titleStr)

end