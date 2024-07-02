% Parameter file for the FSE models
% Initialize parameters automatically when the model gets loaded:

load('FSE_M85.txt');  % load motor characteristics to workspace
load('FSE_M100.txt');


i_Diff = 6;                   % mechanical transmission

TC_On = 1;                    % enable traction control

TC_Derivative = 0.4;          % threshold traction control mode 1
TC_Relay_On = 0.07;           % switching on traction control mode 2
TC_Relay_Off = 0.065;         % switching off traction control mode 2

TC_Incr_Gas = 0.005;          % increasing electronic gas pedal (mode 0)
TC_Decr_Gas = 0.00005;        % decreasing electronic gas pedal (mode 2)

% ------------------------------
% Torque Vectoring Model
% ------------------------------
TV_On = 1;                    % enable torque vectoring
TV_Override = 0.1;              % increase this parameter to increase the torque vectoring effect; can lead to instability

% ------------------------------
% Recuperation Model
% ------------------------------
eta_Gen = -1;                 % efficency of the electric motor in generator mode; always has to negative for the model to work properly
Gen_Loss = 0.15;              % energy loss due to building up a strator torque for braking
Brake_Pedal_Travel = 0.9;     % define after how much brake pedal a hydraulic brake has to be actuated, can be found in the FSAE rules
Brake_Mech_Rec_Ratio = 0.9;   % define how much recuperational braking is still active after the hydralic brake travel is actuated

% ------------------------------
% Block Battery
% ------------------------------
SOC_Init =  50;  % Initial SOC [%]
C_Batt   = 10;   % Battery capacity [kWh]
U0_Batt = 100;   % No-load voltage [V]
R0_Batt = 0.1;   % Resistance [Ohm]
C1_Batt = 15;    % Capacity [F]
R1_Batt = 0.1;   % Resistance [Ohm]
C2_Batt = 15;    % Capacity [F]
R2_Batt = 0.1;   % Resistance [Ohm]

display('loading m-file complete')