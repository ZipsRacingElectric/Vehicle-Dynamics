%{
%% Vehicle Object
Vehicle object to hold vehicle design parameters
Ben Model
9/7/24

Adapted to Simulink models 12/21/24 - Brian Glen

%% Overview:
This object contains all simulation parameters of the vehicle
Use this class to instanciate a VEHICLE object that you will
assign relevant parameters. Use the VEHICLE object to recall
these parameters easily in your simulation

%% Usage:
create new one with:
ZR25 = vehicle('parameters file path and name')

Access variables like:
mass = ZR25.mass

%% Notes:
- This object is instantiated using a spreadsheet containing the static
parameters of the vehicle. When MATLAB or Simulink simulations are run, the
car object is created and vehicle data is pulled from the spreadsheet.
Update the spreadsheet to update vehicle parameters.

- The spreadsheet must have a tab named "parameters". The first column must
have the parameter name as written in the setParameters() function below.
The second column must have the parameter value in the proper unit listed
below.

%% TODO:
- update aero force calculation and parameters for a more complete picture
- load transfer due to aero forces
- setParameters() should be automated

%% Units:
- distance: m
- weight: kg
- force: N
- angles: degrees (only used for static suspension settings. Models must
convert this to radians)
- inertia: kg * m^2
- stiffness: N/m (linear), Nm/deg (angular)
- velocity: m/s
- area: m^2
- frequency: Hz
- acceleration: m/s^2
- density: kg/m^3
%}

classdef vehicle < handle

    properties (Access=private)
        % no touchy!
        parameters
    end

    properties (Access=public)

        %% Static parameters - pulled from excel sheet when car object is created
        % Dimensions
        wheelbase                                                           % m
        track_width_front                                                   % m
        track_width_rear                                                    % m

        % Mass
        vehicle_mass                                                        % kg
        driver_mass                                                         % kg
        corner_mass_front                                                   % kg, unsprun corner mass. Include half of the control arm masses.
        corner_mass_rear                                                    % kg, unsprun corner mass. Include half of the control arm masses.
        front_mass_distribution                                             % percentage on front axle
        cg_height                                                           % m
        yaw_polar_inertia                                                   % kg * m^2, about the yaw (vertical) axis  of the C.G.

        % Tires
        tire_loaded_radius                                                  % m
        gear_ratio                                                          % (# input rotations / # output rotations)
        tire_mu                                                             % the hot tire friction from TTC data                           
        tire_mu_correction_factor                                           % this is a correction factor to account for reduced tire-road friction compared to TTC data. Typically 2/3 but should be tuned if raining or cold outside
        tire_stiffness                                                      % N/m
        tire_width                                                          % m

        % Gearbox / Drivetrain
        Jm
        Js1
        J1
        J2
        Js2
        Jw
        Dm
        D1
        D2
        Dw
        Ks1
        Ks2
        
        % Kinematics
        static_toe_front                                                    % degrees (per wheel), + is toe out
        static_toe_rear                                                     % degrees (per wheel), + is toe out
        static_camber_front                                                 % degrees, - is leaning torwards car
        static_camber_rear                                                  % degrees, - is leaning torwards car
        steering_ratio                                                      % ratio, steering wheel angle / ackerman steering angle (aka average of L and R angles)
        ackermann_percentage                                                % percentage
        steering_arm_length                                                 % m, perpendicular length between tire rod point and kingpin axis
        steering_rack_length                                                % m, eye to eye length of steering rack
        tie_rod_length_front                                               % m, perpendicular length between tire rod point and kingpin axis
        steering_rack_to_axis_distance                                      % m, distance between kingpin axis and steering rack, parallel to the longitudinal plane of the vehicle
        steering_pinion_radius                                              % m, radius of the steering rack pinion gear (reference for gear ratio calculation)
        roll_center_front                                                   % m, height of front roll center at static ride height
        roll_center_rear                                                    % m, height of rear roll center at static ride height
        
        % Aerodyamics
        frontal_area                                                        % m^2
        Cd                                                                  % unitless
        Cl                                                                  % unitless. Certain models may require it to be negative or positive based on implementation
        center_of_pressure_distribution                                     % ratio 0(at rear axle) to 1(at front axle)
        velocity_skidpad                                                    % velocity of skidpad for aero measurement
        cla_at_skidpad                                                      % unitless, + is downforce, ClA at skidpad
        cop_at_skidpad                                                      % ratio 0(at rear axle) to 1(at front axle), Cop at skidpad
        velocity_max                                                        % maximum velocity for aero measurement
        cla_at_max_velocity                                                 % unitless, + is downforce,cla at max velocity
        cop_at_max_velocity                                                 % ratio 0(at rear axle) to 1(at front axle), CoP at max velocity

        % Springs and Dampers
        damper_travel                                                       % m, maximum travel of the damper
        spring_rate_front                                                   % N/m, spring rate at the damper
        spring_rate_rear                                                    % N/m, spring rate at the damper
        bar_spring_rate_front                                               % N/m, Spring rate of front roll bar
        bar_spring_rate_rear                                                % N/m, Spring rate of rear roll bar
        motion_ratio_front                                                  % unitless, Damper / wheel (assumes we use coil-overs)
        motion_ratio_rear                                                   % unitless, Damper / wheel (assumes we use coil-overs)
        bar_motion_ratio_front                                              % unitless, Roll bar / wheel (assumes we use coil-overs)
        bar_motion_ratio_rear                                               % unitless, Roll bar / wheel (assumes we use coil-overs)
        ride_frequency_front                                                % Hz, target front ride frequency (compare to calculated)
        ride_frequency_rear                                                 % Hz, target rear ride frequency (compare to calculated)

        % Compliance
        toe_deflection_rear                                                 % deg per 1kN, per wheel, linear toe deflection from Fy forces, from experimental testing

        % Brakes
        piston_radius_front                                                 % m
        piston_radius_rear                                                  % m
        num_pistons_front                                                   % unitless
        num_pistons_rear                                                    % unitless
        pad_friction_front                                                  % unitless, dynamic and static friction are assumed to be the same
        pad_friction_rear                                                   % unitless, dynamic and static friction are assumed to be the same
        max_pedal_force                                                     % N
        disc_radius_front                                                   % m
        disc_radius_rear                                                    % m
        pad_height_front                                                    % m
        pad_height_rear                                                     % m
        mc_diameter_front                                                   % m
        mc_diameter_rear                                                    % m
        balance_bar_ratio_front                                             % 0 to 1
        brake_pedal_motion_ratio                                            % unitless

        %% Calculated parameters - call update method to calculate
        % Dimensions
        a                                                                   % m, distance from front axle to CG
        b                                                                   % m, distance from rear axle to CG
        b_sprung                                                            % m, distance from rear axle to sprung CG
        a_sprung                                                            % m, distance from front axle to sprung CG
        cg_sprung                                                           % m, height of sprung CG above ground
        average_track_width                                                 % m, average track width (front and rear)
        cg_track_ratio                                                      % unitless, ratio of CG height to mean track wid

        % Mass
        mass_total                                                          % kg, total vehicle mass
        sprung_mass_total                                                   % kg, total sprung mass
        sprung_mass_front                                                   % kg, front axle sprung mass
        sprung_mass_rear                                                    % kg, rear axle sprung mass
        unsprung_mass_total                                                 % kg, total unsprung mass
        unsprung_mass_front                                                 % kg, front axle unsprung mass
        unsprung_mass_rear                                                  % kg, rear axle unsprung mass
        average_corner_mass                                                 % kg, average mass at each corner
        
        front_mass                                                        % kg, total front axle weight
        rear_mass                                                         % kg, total rear axle weight
        front_corner_mass                                                 % kg, average front corner weight
        rear_corner_mass                                                  % kg, average rear corner weight

        % Springs and Dampers
        wheel_rate_front                                                    % N/m, effective front wheel rate
        wheel_rate_rear                                                     % N/m, effective rear wheel rate
        wheel_rate_from_bar_front                                           % N/m, contribution of front anti-roll bar to wheel rate
        wheel_rate_from_bar_rear                                            % N/m, contribution of rear anti-roll bar to wheel rate
        spring_rate_total_front                                             % N/m, combined front spring rate
        spring_rate_total_rear                                              % N/m, combined rear spring rate
        roll_wheel_rate_front                                               % N/m, total front roll wheel rate
        roll_wheel_rate_rear                                                % N/m, total rear roll wheel rate
        total_roll_rate_front                                               % N/m, total front roll rate (including tires)
        total_roll_rate_rear                                                % N/m, total rear roll rate (including tires)

        % Aerodynamics
        downforce_at_skidpad_front                                          % N, aerodynamic downforce on front axle at skidpad velocity
        downforce_at_skidpad_rear                                           % N, aerodynamic downforce on rear axle at skidpad velocity
        downforce_at_max_velocity_front                                     % N, aerodynamic downforce on front axle at max velocity
        downforce_at_max_velocity_rear                                      % N, aerodynamic downforce on rear axle at max velocity

        % Lateral Load Transfer Distribution
        height_nra_to_sm_center                                             % m, height of neutral roll axis relative to sprung mass CG
        kf                                                                  % Nm/deg, front roll stiffness
        kr                                                                  % Nm/deg, rear roll stiffness
        kf_prime                                                            % Nm/deg, adjusted front roll stiffness
        kr_prime                                                            % Nm/deg, adjusted rear roll stiffness
        lltd_front                                                          % percent, lateral load transfer distribution on front axle
        lltd_rear                                                           % percent, lateral load transfer distribution on rear axle
        tlltd                                                               % percent, total lateral load transfer distribution

        % Longitudinal Load Transfer
        long_load_transfer                                                  % N, longitudinal load transfer force

        % Gearbox Dynamics
        Jeq
        Deq
        Keq

        %% Constants
        % Environment
        air_temp = 20;                                                      % degrees C
        air_density = 1.225;                                                % kg / m^3, at 20 degrees C
        g = 9.81;                                                           % m/s^2, acceleration of gravity
    end

    methods

        %{

        Creates a vehicle object holding vehicle parameters.
        Pass it a spreadsheet containing vehicle data (see above)

        %}
        function obj = vehicle(parametersFile)
            % UNTITLED Construct an instance of this class
            % Constructor function - pass parameters file
            parameterTable = readtable(parametersFile,'Sheet','parameters');
            
            % stores unchanged parameters from Excel file for reset
            % private variable, no touchy!
            obj.parameters = containers.Map(parameterTable.Var1, parameterTable.Var2);
            
            % set non-calculated variables in class to values from
            % parameters container
            obj.setParameters();            

            % evaluate calculated parameters
            obj.update(); 

            disp('Vehicle object created successfully!');

        end

        %{
        
        Reads raw values from the excel obect and assigns values to the vehicle object parameters
        
        %}
        function setParameters(obj)

            %% Set Static Parameters
            % Dimensions
            obj.wheelbase = obj.parameters('wheelbase');                  
            obj.track_width_front = obj.parameters('track_width_front');   
            obj.track_width_rear = obj.parameters('track_width_rear');      
        
            % Mass
            obj.vehicle_mass = obj.parameters('vehicle_mass');          
            obj.driver_mass = obj.parameters('driver_mass');                   
            obj.corner_mass_front = obj.parameters('corner_mass_front');          
            obj.corner_mass_rear = obj.parameters('corner_mass_rear');            
            obj.front_mass_distribution = obj.parameters('front_mass_distribution'); 
            obj.cg_height = obj.parameters('cg_height');                    
            obj.yaw_polar_inertia = obj.parameters('yaw_polar_inertia');    
        
            % Tires
            obj.tire_loaded_radius = obj.parameters('tire_loaded_radius');     
            obj.gear_ratio = obj.parameters('gear_ratio');                 
            obj.tire_mu = obj.parameters('tire_mu');     
            obj.tire_mu_correction_factor = obj.parameters('tire_mu_correction_factor')
            obj.tire_stiffness = obj.parameters('tire_stiffness');        
            obj.tire_width = obj.parameters('tire_width');

            % Gearbox / Drivetrain
            obj.Jm = obj.parameters('Jm');
            obj.Js1 = obj.parameters('Js1');
            obj.J1 = obj.parameters('J1');
            obj.J2 = obj.parameters('J2');
            obj.Js2 = obj.parameters('Js2');
            obj.Jw = obj.parameters('Jw');
            obj.Dm = obj.parameters('Dm');
            obj.D1 = obj.parameters('D1');
            obj.D2 = obj.parameters('D2');
            obj.Dw = obj.parameters('Dw');
            obj.Ks1 = obj.parameters('Ks1');
            obj.Ks2 = obj.parameters('Ks2');
        
            % Kinematics
            obj.static_toe_front = obj.parameters('static_toe_front');         
            obj.static_toe_rear = obj.parameters('static_toe_rear');             
            obj.static_camber_front = obj.parameters('static_camber_front');      
            obj.static_camber_rear = obj.parameters('static_camber_rear');        
            obj.steering_ratio = obj.parameters('steering_ratio');                
            obj.ackermann_percentage = obj.parameters('ackermann_percentage');      
            obj.steering_arm_length = obj.parameters('steering_arm_length');      
            obj.steering_rack_length = obj.parameters('steering_rack_length');     
            obj.tie_rod_length_front = obj.parameters('tie_rod_length_front');  
            obj.steering_rack_to_axis_distance = obj.parameters('steering_rack_to_axis_distance');
            obj.steering_pinion_radius = obj.parameters('steering_pinion_radius');  
            obj.roll_center_front = obj.parameters('roll_center_front');        
            obj.roll_center_rear = obj.parameters('roll_center_rear');        
        
            % Aerodynamics
            obj.frontal_area = obj.parameters('frontal_area');                 
            obj.Cd = obj.parameters('Cd');                                    
            obj.Cl = obj.parameters('Cl');                                     
            obj.center_of_pressure_distribution = obj.parameters('center_of_pressure_distribution'); 
            obj.velocity_skidpad = obj.parameters('velocity_skidpad');            
            obj.cla_at_skidpad = obj.parameters('cla_at_skidpad');               
            obj.cop_at_skidpad = obj.parameters('cop_at_skidpad');                  
            obj.velocity_max = obj.parameters('velocity_max');                     
            obj.cla_at_max_velocity = obj.parameters('cla_at_max_velocity');       
            obj.cop_at_max_velocity = obj.parameters('cop_at_max_velocity');        
        
            % Springs and Dampers
            obj.damper_travel = obj.parameters('damper_travel');             
            obj.spring_rate_front = obj.parameters('spring_rate_front');          
            obj.spring_rate_rear = obj.parameters('spring_rate_rear');            
            obj.bar_spring_rate_front = obj.parameters('bar_spring_rate_front');    
            obj.bar_spring_rate_rear = obj.parameters('bar_spring_rate_rear');     
            obj.motion_ratio_front = obj.parameters('motion_ratio_front');       
            obj.motion_ratio_rear = obj.parameters('motion_ratio_rear');            
            obj.bar_motion_ratio_front = obj.parameters('bar_motion_ratio_front');  
            obj.bar_motion_ratio_rear = obj.parameters('bar_motion_ratio_rear');    
            obj.ride_frequency_front = obj.parameters('ride_frequency_front');    
            obj.ride_frequency_rear = obj.parameters('ride_frequency_rear');

            % Compliance
            obj.toe_deflection_rear = obj.parameters('toe_deflection_rear');

            % Brakes
            obj.piston_radius_front = obj.parameters('piston_radius_front');
            obj.piston_radius_rear = obj.parameters('piston_radius_rear');
            obj.num_pistons_front = obj.parameters('num_pistons_front');
            obj.num_pistons_rear = obj.parameters('num_pistons_rear');
            obj.pad_friction_front = obj.parameters('pad_friction_front');
            obj.pad_friction_rear = obj.parameters('pad_friction_rear');
            obj.max_pedal_force = obj.parameters('max_pedal_force');
            obj.disc_radius_front = obj.parameters('disc_radius_front');
            obj.disc_radius_rear = obj.parameters('disc_radius_rear');
            obj.pad_height_front = obj.parameters('pad_height_front');
            obj.pad_height_rear = obj.parameters('pad_height_rear');
            obj.mc_diameter_front = obj.parameters('mc_diameter_front');
            obj.mc_diameter_rear = obj.parameters('mc_diameter_rear');
            obj.balance_bar_ratio_front = obj.parameters('balance_bar_ratio_front');
            obj.brake_pedal_motion_ratio = obj.parameters('brake_pedal_motion_ratio');
        end
        
        %{
        
        Updates the calculated parameters from the static parameters. This
        needs to be called everytime the static parameters change.
        
        %}
        function update(obj)

            %% Evaluate calculated parameters
            % Mass
            obj.mass_total = obj.vehicle_mass + obj.driver_mass;
            obj.unsprung_mass_front = 2 * obj.corner_mass_front;
            obj.unsprung_mass_rear = 2 * obj.corner_mass_rear;
            obj.unsprung_mass_total = obj.unsprung_mass_front + obj.unsprung_mass_rear;
            obj.sprung_mass_total = obj.mass_total - obj.unsprung_mass_total;
            obj.sprung_mass_front = (obj.mass_total * obj.front_mass_distribution) - obj.unsprung_mass_front;
            obj.sprung_mass_rear = obj.sprung_mass_total - obj.sprung_mass_front;
            obj.average_corner_mass = obj.unsprung_mass_total / 4;
            obj.front_mass = obj.mass_total * obj.front_mass_distribution;
            obj.rear_mass = obj.mass_total * (1 - obj.front_mass_distribution);
            obj.front_corner_mass = obj.front_mass / 2;                             % total weight (in kg) on a front wheel
            obj.rear_corner_mass = obj.rear_mass / 2;                               % total weight (in kg) on a rear wheel
            
            % Dimensions
            obj.a = obj.wheelbase * (1 - obj.front_mass_distribution);
            obj.b = obj.wheelbase * (obj.front_mass_distribution);          
            obj.b_sprung = (obj.mass_total * obj.b - obj.unsprung_mass_front * obj.wheelbase) / obj.sprung_mass_total;
            obj.a_sprung = obj.wheelbase - obj.b_sprung;
            obj.cg_sprung = ((obj.mass_total * obj.cg_height) - (obj.unsprung_mass_front * obj.tire_loaded_radius) - (obj.unsprung_mass_rear * obj.tire_loaded_radius)) / obj.sprung_mass_total;
            obj.average_track_width = (obj.track_width_front + obj.track_width_rear) / 2;
            obj.cg_track_ratio = obj.cg_height / obj.average_track_width;
            
            % Spring and Dampers
            obj.wheel_rate_front = obj.spring_rate_front * obj.motion_ratio_front^2;
            obj.wheel_rate_rear = obj.spring_rate_rear * obj.motion_ratio_rear^2;
            obj.wheel_rate_from_bar_front = obj.bar_spring_rate_front * obj.bar_motion_ratio_front^2;
            obj.wheel_rate_from_bar_rear = obj.bar_spring_rate_rear * obj.bar_motion_ratio_rear^2;
            obj.spring_rate_total_front = (obj.wheel_rate_front * obj.tire_stiffness) / (obj.wheel_rate_front + obj.tire_stiffness);
            obj.spring_rate_total_rear = (obj.wheel_rate_rear * obj.tire_stiffness) / (obj.wheel_rate_rear + obj.tire_stiffness);
            obj.roll_wheel_rate_front = obj.wheel_rate_front + obj.wheel_rate_from_bar_front;
            obj.roll_wheel_rate_rear = obj.wheel_rate_rear + obj.wheel_rate_from_bar_rear;
            obj.total_roll_rate_front = (obj.roll_wheel_rate_front * obj.tire_stiffness) / (obj.roll_wheel_rate_front + obj.tire_stiffness);
            obj.total_roll_rate_rear = (obj.roll_wheel_rate_rear * obj.tire_stiffness) / (obj.roll_wheel_rate_rear + obj.tire_stiffness);

            % Aerodynamics
            obj.downforce_at_skidpad_front = (obj.cla_at_skidpad * obj.air_density * 0.5 * obj.velocity_skidpad^2) * obj.center_of_pressure_distribution;
            obj.downforce_at_skidpad_rear = (obj.cla_at_skidpad * obj.air_density * 0.5 * obj.velocity_skidpad^2) * (1 - obj.center_of_pressure_distribution);
            obj.downforce_at_max_velocity_front = (obj.cla_at_max_velocity * obj.air_density * 0.5 * obj.velocity_max^2) * obj.cop_at_max_velocity;
            obj.downforce_at_max_velocity_rear = (obj.cla_at_max_velocity * obj.air_density * 0.5 * obj.velocity_max^2) * (1 - obj.cop_at_max_velocity);

            % Lateral Load Transfer Calculations
            obj.height_nra_to_sm_center = abs((obj.roll_center_rear - obj.roll_center_front) * (obj.a - obj.a_sprung) - (-obj.b - obj.a) * obj.cg_sprung + (-obj.b) * obj.roll_center_front - obj.roll_center_rear * obj.a) / sqrt((obj.roll_center_front - obj.roll_center_rear)^2 + (-obj.b - obj.a)^2);
            obj.kf = 0.5 * obj.total_roll_rate_front * obj.track_width_front^2;
            obj.kr = 0.5 * obj.total_roll_rate_rear * obj.track_width_rear^2;
            obj.kf_prime = obj.kf - (obj.wheelbase - obj.a_sprung) * obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center / obj.wheelbase;
            obj.kr_prime = obj.kr - obj.a_sprung * obj.sprung_mass_total * obj.height_nra_to_sm_center / obj.wheelbase;
            [obj.lltd_front, obj.lltd_rear] = obj.calculate_lateral_load_transfer_front();
            obj.tlltd = obj.lltd_front / (obj.lltd_front + obj.lltd_rear);

            % Longitudinal Load Transfer Calculations
            obj.long_load_transfer = obj.cg_height / obj.wheelbase * obj.mass_total * obj.g;

            % Gearbox Dynamics
            obj.Jeq = (obj.J1 + obj.Js1 + obj.Jm) + (obj.J2 + obj.Js2 + obj.Jw) * (1 / obj.gear_ratio)^2;
            obj.Deq = (obj.D1 + obj.Dm) + (obj.D2 + obj.Dw) * (1 / obj.gear_ratio)^2;
            obj.Keq = obj.Ks1 + obj.Ks2 * (1 / obj.gear_ratio)^2;
        end

        %{
        
        Sets all parameters back to default from excel file

        %}
        function reset(obj)
            obj.setParameters()
            obj.update()
        end
        
        %{
        
        Creates Simulink parameters in the workspace for each public parameter
    
        %}
        function create_simulink_parameters(obj)

            % Get all public properties of the object
            parameter_names = properties(obj);
            
            if isempty(parameter_names)
                disp('No public properties found in the object.');
                return;
            end

            % Iterate through each parameter
            for i = 1:length(parameter_names)
                param_name = parameter_names{i};
                param_value = obj.(param_name);

                % Skip empty properties
                if isempty(param_value)
                    warning(['Parameter ', param_name, ' is empty. Skipping creation of simulink parameter.']);
                    continue;
                end
                
                % Create the simulink parameter in the base workspace
                assignin('base', param_name, Simulink.Parameter(param_value));
                fprintf('Created the parameter %s = %d\n', param_name, param_value);
            end
            
            disp('Simulink parameters created in the workspace.');
        end
        
        %% Helper Functions for vehicle analysis

        %{

        Calculates lateral load transfers

        %}
        function [lltd_front, lltd_rear] = calculate_lateral_load_transfer_front(obj)
            % Calculate the lateral load transfer distribution for the front axle
            front_roll_center = obj.roll_center_front;
            rear_roll_center = obj.roll_center_rear;
            height_sprung_mass_cg = obj.cg_sprung;

            lltd_front = (obj.sprung_mass_total * obj.g / obj.track_width_front) * ...
                ((obj.height_nra_to_sm_center * obj.kf_prime) / (obj.kf + obj.kr - obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center) + ...
                (obj.wheelbase - obj.a_sprung) / obj.wheelbase * front_roll_center) + ...
                obj.unsprung_mass_front * obj.g / obj.track_width_front * height_sprung_mass_cg;

            lltd_rear = (obj.sprung_mass_total * obj.g / obj.track_width_rear) * ...
                ((obj.height_nra_to_sm_center * obj.kr_prime) / (obj.kf + obj.kr - obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center) + ...
                (obj.a_sprung) / obj.wheelbase * rear_roll_center) + ...
                obj.unsprung_mass_rear * obj.g / obj.track_width_rear * height_sprung_mass_cg;
        end
        
        %{
        
        Calculates aero forces in [N] at a given velocity
        
        %}
        function [downforce, dragforce] = calculate_aero_forces(obj, velocity)
            % Calculate the drag force given a velocity
            dragforce = 0.5 * obj.air_density * velocity^2 * obj.Cd * obj.frontal_area;
            % Calculate the aerodynamic downforce given a velocity
            downforce = 0.5 * obj.air_density * velocity^2 * obj.Cl * obj.frontal_area;
        end

    end
end