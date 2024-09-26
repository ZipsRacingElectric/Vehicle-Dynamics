% Vehicle parameters for ZR25 design.
% 9/26/24

% Dimensions
wheelbase = 1.530;                              % m
front_track_width = 1.230;                      % m
rear_track_width = 1.2;                         % m

% Mass
vehicle_mass = 205;                             % kg
driver_mass = 72;                               % kg
front_mass_distribution = 0.48;                 % percentage on front axle
cg_height = 0.255;                              % m
a = wheelbase * (1 - front_mass_distribution);  % m
b = wheelbase * (front_mass_distribution);      % m
yaw_polar_inertia = 80;                         % kg * m^2

% Tires
tire_loaded_radius = 0.25;                      % m
gear_ratio = 13;                                % (# input rotations / # output rotations)
tire_mu = 2;                                    % use only if tire model not availible                           
tire_stiffness = 76000;                         % N/m
tire_width = 0.19;                              % m

% Aerodyamics
frontal_area = 1.5;                             % m^2
CdA = 0.2;                                      % unitless
ClA = -3;                                       % unitless. Certain models may require it to be negative or positive based on implementation
center_of_pressure = 0.45;                      % percentage from front axle

% Kinematics
static_toe_front = 0.5;                         % degrees (per wheel), + is toe out
static_toe_rear = 0;                            % degrees (per wheel), + is toe out
static_camber_front = -2.5;                     % degrees, - is leaning torwards car
static_camber_rear = -2.5;                      % degrees, - is leaning torwards car

% Springs and Dampers
spring_rate_front = 35000;                      % N/m
spring_rate_rear = 40000;                       % N/m
motion_ratio_front = 1;                         % Damper / Wheel (assumes we use coil-overs)
motion_ratio_rear = 1;                          % Damper / Wheel (assumes we use coil-overs)
