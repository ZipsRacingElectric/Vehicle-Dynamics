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
- create a function to build Simulink parameters in the workspace

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
        tire_mu                                                             % use only if tire model not availible                           
        tire_stiffness                                                      % N/m
        tire_width                                                          % m
        
        % Kinematics
        static_toe_front                                                    % degrees (per wheel), + is toe out
        static_toe_rear                                                     % degrees (per wheel), + is toe out
        static_camber_front                                                 % degrees, - is leaning torwards car
        static_camber_rear                                                  % degrees, - is leaning torwards car
        steering_ratio                                                      % ratio, steering wheel angle / ackerman steering angle (aka average of L and R angles)
        ackermann_percentage                                                % percentage
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

        % these are sorta references
        ride_frequency_front                                                % Hz, target front ride frequency (compare to calculated)
        ride_frequency_rear                                                 % Hz, target rear ride frequency (compare to calculated)

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
        
        front_weight                                                        % kg, total front axle weight
        rear_weight                                                         % kg, total rear axle weight
        front_corner_weight                                                 % kg, average front corner weight
        rear_corner_weight                                                  % kg, average rear corner weight

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

        %% Constants
        % Environment
        air_density = 1.225;                                                % kg / m^3, at 20 degrees C
        air_density_30C = 1.164;                                            % kg / m^3, at 30 degrees C
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
            obj.wheelbase = obj.parameters('wheelbase');                              % m
            obj.track_width_front = obj.parameters('track_width_front');              % m
            obj.track_width_rear = obj.parameters('track_width_rear');                % m
        
            % Mass
            obj.vehicle_mass = obj.parameters('vehicle_mass');                        % kg
            obj.driver_mass = obj.parameters('driver_mass');                          % kg
            obj.corner_mass_front = obj.parameters('corner_mass_front');              % kg
            obj.corner_mass_rear = obj.parameters('corner_mass_rear');                % kg
            obj.front_mass_distribution = obj.parameters('front_mass_distribution');  % percentage
            obj.cg_height = obj.parameters('cg_height');                              % m
            obj.yaw_polar_inertia = obj.parameters('yaw_polar_inertia');              % kg·m²
        
            % Tires
            obj.tire_loaded_radius = obj.parameters('tire_loaded_radius');            % m
            obj.gear_ratio = obj.parameters('gear_ratio');                            % unitless
            obj.tire_mu = obj.parameters('tire_mu');                                  % unitless
            obj.tire_stiffness = obj.parameters('tire_stiffness');                    % N/m
            obj.tire_width = obj.parameters('tire_width');                            % m
        
            % Kinematics
            obj.static_toe_front = obj.parameters('static_toe_front');                % degrees
            obj.static_toe_rear = obj.parameters('static_toe_rear');                  % degrees
            obj.static_camber_front = obj.parameters('static_camber_front');          % degrees
            obj.static_camber_rear = obj.parameters('static_camber_rear');            % degrees
            obj.steering_ratio = obj.parameters('steering_ratio');                    % unitless
            obj.ackermann_percentage = obj.parameters('ackermann_percentage');          % percentage
            obj.roll_center_front = obj.parameters('roll_center_front');              % m
            obj.roll_center_rear = obj.parameters('roll_center_rear');                % m
        
            % Aerodynamics
            obj.frontal_area = obj.parameters('frontal_area');                        % m²
            obj.Cd = obj.parameters('Cd');                                            % unitless
            obj.Cl = obj.parameters('Cl');                                            % unitless
            obj.center_of_pressure_distribution = obj.parameters('center_of_pressure_distribution'); % unitless
            obj.velocity_skidpad = obj.parameters('velocity_skidpad');                % m/s
            obj.cla_at_skidpad = obj.parameters('cla_at_skidpad');                    % unitless
            obj.cop_at_skidpad = obj.parameters('cop_at_skidpad');                    % unitless
            obj.velocity_max = obj.parameters('velocity_max');                        % m/s
            obj.cla_at_max_velocity = obj.parameters('cla_at_max_velocity');          % unitless
            obj.cop_at_max_velocity = obj.parameters('cop_at_max_velocity');          % unitless
        
            % Springs and Dampers
            obj.damper_travel = obj.parameters('damper_travel');                      % m
            obj.spring_rate_front = obj.parameters('spring_rate_front');              % N/m
            obj.spring_rate_rear = obj.parameters('spring_rate_rear');                % N/m
            obj.bar_spring_rate_front = obj.parameters('bar_spring_rate_front');      % N/m
            obj.bar_spring_rate_rear = obj.parameters('bar_spring_rate_rear');        % N/m
            obj.motion_ratio_front = obj.parameters('motion_ratio_front');            % unitless
            obj.motion_ratio_rear = obj.parameters('motion_ratio_rear');              % unitless
            obj.bar_motion_ratio_front = obj.parameters('bar_motion_ratio_front');    % unitless
            obj.bar_motion_ratio_rear = obj.parameters('bar_motion_ratio_rear');      % unitless
        
            % Ride Frequencies
            obj.ride_frequency_front = obj.parameters('ride_frequency_front');        % Hz
            obj.ride_frequency_rear = obj.parameters('ride_frequency_rear');          % Hz
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
            obj.front_weight = obj.mass_total * obj.front_mass_distribution;
            obj.rear_weight = obj.mass_total * (1 - obj.front_mass_distribution);
            obj.front_corner_weight = obj.front_weight / 2;
            obj.rear_corner_weight = obj.rear_weight / 2;
            
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
            obj.downforce_at_skidpad_front = (obj.cla_at_skidpad * obj.air_density * 0.5 * obj.velocity_skidpad^2) * obj.center_of_pressure_distribution / obj.g;
            obj.downforce_at_skidpad_rear = (obj.cla_at_skidpad * obj.air_density * 0.5 * obj.velocity_skidpad^2) * (1 - obj.center_of_pressure_distribution) / obj.g;
            obj.downforce_at_max_velocity_front = (obj.cla_at_max_velocity * obj.air_density * 0.5 * obj.velocity_max^2) * obj.cop_at_max_velocity / obj.g;
            obj.downforce_at_max_velocity_rear = (obj.cla_at_max_velocity * obj.air_density * 0.5 * obj.velocity_max^2) * (1 - obj.cop_at_max_velocity) / obj.g;

            % Lateral Load Transfer Calculations
            obj.height_nra_to_sm_center = abs((obj.roll_center_rear - obj.roll_center_front) * (obj.a - obj.a_sprung) - (-obj.b - obj.a) * obj.cg_sprung + (-obj.b) * obj.roll_center_front - obj.roll_center_rear * obj.a) / sqrt((obj.roll_center_front - obj.roll_center_rear)^2 + (-obj.b - obj.a)^2);
            obj.kf = 0.5 * obj.total_roll_rate_front * obj.track_width_front^2;
            obj.kr = 0.5 * obj.total_roll_rate_rear * obj.track_width_rear^2;
            obj.kf_prime = obj.kf - (obj.wheelbase - obj.a_sprung) * obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center / obj.wheelbase;
            obj.kr_prime = obj.kr - obj.a_sprung * obj.sprung_mass_total * obj.height_nra_to_sm_center / obj.wheelbase;
            obj.lltd_front = obj.calculate_lateral_load_transfer_front();
            obj.lltd_rear = obj.calculate_lateral_load_transfer_rear();
            obj.tlltd = obj.lltd_front / (obj.lltd_front + obj.lltd_rear);

            % Longitudinal Load Transfer Calculations
            obj.long_load_transfer = obj.cg_height / obj.wheelbase * obj.mass_total * obj.g;
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

                % Create a Simulink.Parameter object for each property
                
                % Create the simulink parameter in the base workspace
                assignin('base', param_name, Simulink.Parameter(param_value));
                fprintf('Created the parameter %s = %d\n', param_name, param_value);
            end
            
            disp('Simulink parameters created in the workspace.');
        end
        
        %% Helper Functions for vehicle analysis

        function lltd_front = calculate_lateral_load_transfer_front(obj)
            % Calculate the lateral load transfer distribution for the front axle
            front_roll_center = obj.roll_center_front;
            height_sprung_mass_cg = obj.cg_sprung;

            lltd_front = (obj.sprung_mass_total * obj.g / obj.track_width_front) * ...
                ((obj.height_nra_to_sm_center * obj.kf_prime) / (obj.kf + obj.kr - obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center) + ...
                (obj.wheelbase - obj.a_sprung) / obj.wheelbase * front_roll_center) + ...
                obj.unsprung_mass_front * obj.g / obj.track_width_front * height_sprung_mass_cg;
        end
        
        function lltd_rear = calculate_lateral_load_transfer_rear(obj)
            % Calculate the lateral load transfer distribution for the rear axle
            rear_roll_center = obj.roll_center_rear;
            height_sprung_mass_cg = obj.cg_sprung;

            lltd_rear = (obj.sprung_mass_total * obj.g / obj.track_width_rear) * ...
                ((obj.height_nra_to_sm_center * obj.kr_prime) / (obj.kf + obj.kr - obj.sprung_mass_total * obj.g * obj.height_nra_to_sm_center) + ...
                (obj.a_sprung) / obj.wheelbase * rear_roll_center) + ...
                obj.unsprung_mass_rear * obj.g / obj.track_width_rear * height_sprung_mass_cg;
        end
        
        function drag_force = calculate_drag_force(obj, velocity)
            % Calculate the drag force given a velocity
            drag_force = 0.5 * obj.air_density * velocity^2 * obj.Cd * obj.frontal_area;
        end

        function downforce = calculate_downforce(obj, velocity)
            % Calculate the aerodynamic downforce given a velocity
            downforce = 0.5 * obj.air_density * velocity^2 * obj.Cl * obj.frontal_area;
        end

    end
end