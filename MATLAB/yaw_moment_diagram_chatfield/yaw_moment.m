%{

TV Yaw Rate calculator
by Chris Chatfield - 2022-03-06 % Updated 2023-02-12

%% Overview
The yaw moment diagram code from Chatfield's master thesis has been adapted
to use the vehicle object.

%% TODO
- TW_to_SR function is missing from masters thesis, need to write function
- static toe and camber
- better aerodynamic forces
- compliance (toe, camber, but mainly steering compliance)

%}
clear, clc;

% ensure all relevent folders are on the MATLAB path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

%% Load Vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

%% Initial vehicle constants
g = zr25.g;                     % gravitational
Ut = zr25.tire_mu;              % tire coeff of friction (not used)
m = zr25.mass_total;            % Total mass of vehicle
W = m * g;                      % weight of car
L = zr25.wheelbase;             % wheelbase (m)
t_f = zr25.track_width_front;   % vehicle track
t_r = zr25.track_width_rear;
total_width = mean([t_f, t_r]) + zr25.tire_width;             % vehicle total
Hcg = zr25.cg_height;                    % CG height (m)
Fwr = zr25.front_mass_distribution;      % Front weight ratio
Rwr = 1 - Fwr;                  % Rear weight ratio
Rt = zr25.tire_loaded_radius;   % tire radius (m)
Ct = 2*pi*Rt;                   % tire circumference (m)
gear = zr25.gear_ratio;         % gearbox gear reduction ratio
B_rs = 0.52;                    % rear roll stiffness fraction
pr = 0.3;                       % relaxation coefficient for yaw iteration

dFzf_dAxb = Hcg * m / (2*L);
dFzf_dAyb = Hcg * W * (1 - B_rs) / t_f; 
dFzr_dAyb = Hcg * W * B_rs / t_r; 

% motor properties
peak_torque = 21;               % peak torque (N-m)
peak_power = 35E3;              % peak power (W)
peak_TW = peak_torque * gear;   % peak torque at wheels after gearing

%% weight transfer calc
Wf = Fwr * W;                   % Weight of front (N)
Wr = Rwr * W;                   % Weight of rear (N)
a = Wr * L / W;                 % CG dist from front axle (m)
b = Wf * L / W;                 % CG dist from rear axle (m)

% tire contact patch coordinates (Ypos to left)
Y_fl_b = t_f/2;
Y_fr_b = -Y_fl_b;
Y_rl_b = t_r/2;
Y_rr_b = -Y_rl_b;

X_fl_b = a;
X_fr_b = X_fl_b;
X_rl_b = -b;
X_rr_b = X_rl_b;

r = 8.35;                       % radius of turn (m)
Cf = 750;                       % front cornering force (N/deg)
Cr = 250;                       % rear cornering force (N/deg)
Vx_mph = 30;                    % body velocity (mph)
Vx = Vx_mph * 0.44704;          % body velocity (m/s)

% aero loads
Fadf = 0.202*(Vx_mph^2);        % Front axle downforce (N)
Radf = 0.317*(Vx_mph^2);        % Rear axle downforce (N)
DF = Fadf + Radf;

% wheel velocities, infinite radius corner
Vx_fl = Vx; % * 7.74 / 8.35;
Vx_fr = Vx; % * 8.96 / 8.35;
Vx_rl = Vx_fl;
Vx_rr = Vx_fr;

%% vehicle slip and steer angle generation, with TV proportional to steering angle delta
tv_to_delta_deg_ratio = 21; %21
tvd = tv_to_delta_deg_ratio * 180/pi();

betas_deg = -11:1:11; % vehicle slip grid
deltas_deg = -15:1:15; % front steer angle grid, assuming rear delta = 0

betas = betas_deg * pi/180;
deltas = deltas_deg * pi/180;

Cn_mat = zeros(size(betas,2),size(deltas,2));
Ax_mat = zeros(size(betas,2),size(deltas,2));
Ay_mat = zeros(size(betas,2),size(deltas,2));
i_mat = zeros(size(betas,2),size(deltas,2));

b = 1;
for beta = betas
    d = 1;

    % convert corner coordinates to velocity frame
    X_fl_v = X_fl_b*cos(beta) + Y_fl_b*sin(beta);
    X_fr_v = X_fr_b*cos(beta) + Y_fr_b*sin(beta);
    X_rl_v = X_rl_b*cos(beta) + Y_rl_b*sin(beta);
    X_rr_v = X_rr_b*cos(beta) + Y_rr_b*sin(beta);

    Y_fl_v = Y_fl_b*cos(beta) - X_fl_b*sin(beta);
    Y_fr_v = Y_fr_b*cos(beta) - X_fr_b*sin(beta);
    Y_rl_v = Y_rl_b*cos(beta) - X_rl_b*sin(beta);
    Y_rr_v = Y_rr_b*cos(beta) - X_rr_b*sin(beta);

    for delta = deltas

        % initial Fz estimate - no weight transfer
        Fz_fl = (Wf + Fadf) / 2;
        Fz_fr = (Wf + Fadf) / 2;
        Fz_rl = (Wr + Radf) / 2;
        Fz_rr = (Wr + Radf) / 2;

        Ay = 0;
        omega = 0;
        % calculate torque command for steering angle
        TW = min(abs(tvd*delta),peak_TW)*sign(delta); % torque at wheels
        %TWs(d) = TW;

        i = 1; % iteration count

        % compute force balance between load transfer and tire forces
        converge = 0;
        while converge == 0 && i<100
            
            % calculate beta contribution to tire slip angles
            beta_fl = atan((Vx*sin(beta) + omega*X_fl_v) / (Vx*cos(beta) - omega*Y_fl_v));
            beta_fr = atan((Vx*sin(beta) + omega*X_fr_v) / (Vx*cos(beta) - omega*Y_fr_v));
            beta_rl = atan((Vx*sin(beta) + omega*X_rl_v) / (Vx*cos(beta) - omega*Y_rl_v));
            beta_rr = atan((Vx*sin(beta) + omega*X_rr_v) / (Vx*cos(beta) - omega*Y_rr_v));

            % slip angle on each tire
            alpha_fl = beta_fl - delta;
            alpha_fr = beta_fr - delta;
            alpha_rl = beta_rl;
            alpha_rr = beta_rr;

            % slip ratio on each tire based on TW and normal load
            % Note: TW_to_SR function is missing from masters thesis,
            % assuming SR = 0 and no torque vectoring for now
            %{
            SR_fl = TW_to_SR(-TW,Fz_fl,Vx_fl,alpha_fl,'l');
            SR_fr = TW_to_SR( TW,Fz_fr,Vx_fr,alpha_fr,'r');
            SR_rl = TW_to_SR(-TW,Fz_rl,Vx_rl,alpha_rl,'l');
            SR_rr = TW_to_SR( TW,Fz_rr,Vx_rr,alpha_rr,'r');
            %}

            SR_fl = 0;
            SR_fr = 0;
            SR_rl = 0;
            SR_rr = 0;

            
            %{
            if Fz_fl < 1 || Fz_fr < 1 || Fz_rl < 1 || Fz_rr < 1
                SR_fl = 0.1*sign(SR_fl);
                SR_fr = 0.1*sign(SR_fr);
                SR_rl = 0.1*sign(SR_rl);
                SR_rr = 0.1*sign(SR_rr);
            end
            %}

            % apply tire model to get F&M's for each tire
            [Fx_fl,Fy_fl,Mx_fl,My_fl,Mz_fl] = GY_20x7_13_7in_12psi_V43(Fz_fl,Vx_fl,0,alpha_fl,SR_fl,'l');
            [Fx_fr,Fy_fr,Mx_fr,My_fr,Mz_fr] = GY_20x7_13_7in_12psi_V43(Fz_fr,Vx_fr,0,alpha_fr,SR_fr,'r');
            [Fx_rl,Fy_rl,Mx_rl,My_rl,Mz_rl] = GY_20x7_13_7in_12psi_V43(Fz_rl,Vx_rl,0,alpha_rl,SR_rl,'l');
            [Fx_rr,Fy_rr,Mx_rr,My_rr,Mz_rr] = GY_20x7_13_7in_12psi_V43(Fz_rr,Vx_rr,0,alpha_rr,SR_rr,'r');

            % total vehicle forces
            Fx = Fx_fl*cos(delta) - Fy_fl*sin(delta) + Fx_fr*cos(delta) - Fy_fr*sin(delta) + Fx_rl + Fx_rr;
            Fy = Fx_fl*sin(delta) + Fy_fl*cos(delta) + Fx_fr*sin(delta) + Fy_fr*cos(delta) + Fy_rl + Fy_rr;
            Mz = X_fl_v*Fy_fl*cos(delta) + X_fl_v*Fx_fl*sin(delta) - Y_fl_v*Fx_fl*cos(delta) + Y_fl_v*Fy_fl*sin(delta) + X_fr_v*Fx_fr*sin(delta) - Y_fr_v*Fx_fr*cos(delta) + Y_fr_v*Fy_fr*sin(delta) + X_fr_v*Fy_fr*cos(delta) + -Y_rl_v*Fx_rl + X_rl_v*Fy_rl + -Y_rr_v*Fx_rr + X_rr_v*Fy_rr;
            
            Ax = Fx / W;
            Ay_new = Fy / W;
            Cn = Mz / (m*g*L);

            % yaw rate iteration estimation with relaxation (ref Chris Patton eq #70)
            omega_est = Ay_new / Vx;
            omega = omega_est*(1-pr) + omega*pr;

            % calculate weight transfer based on rates
            % max statement added to avoid negative corner loads
            Fz_fl = max(Wf/2 + dFzf_dAxb * Ax + dFzf_dAyb * Ay_new + Fadf/2, 0.1);
            Fz_fr = max(Wf/2 + dFzf_dAxb * Ax - dFzf_dAyb * Ay_new + Fadf/2, 0.1);
            Fz_rl = max(Wr/2 - dFzf_dAxb * Ax + dFzr_dAyb * Ay_new + Radf/2, 0.1);
            Fz_rr = max(Wr/2 - dFzf_dAxb * Ax - dFzr_dAyb * Ay_new + Radf/2, 0.1);

            % logic to add loads back to loaded corners proportionally in case of wheel lifting
            Fz_total = Fz_fl + Fz_fr + Fz_rl + Fz_rr;
            if Fz_total > W + DF
                Fz_fl = Fz_fl * (W+DF) / Fz_total;
                Fz_fr = Fz_fr * (W+DF) / Fz_total;
                Fz_rl = Fz_rl * (W+DF) / Fz_total;
                Fz_rr = Fz_rr * (W+DF) / Fz_total;
                Fz_total = Fz_fl + Fz_fr + Fz_rl + Fz_rr;
            end

            if abs(Ay_new - Ay) < abs(Ay_new*.001) || abs(Ay_new) == 0
                converge = 1;
            else
                Ay = Ay_new;
                i = i+1;
            end

            i_mat(b,d) = i;
            Ax_mat(b,d) = Ax;
            Cn_mat(b,d) = Cn;
            Ay_mat(b,d) = Ay_new;
        end

        d = d+1;
    end
    b = b+1;
end

% plot diagram
figure()
xlabel('Lateral Acceleration Ay')
ylabel('Yaw Moment Cn')
zlabel('Longitudinal Accleration Ax')
hold on

grid on
for d = 1:size(deltas,2)
    l1 = plot3(Ay_mat(:,d),Cn_mat(:,d),Ax_mat(:,d),'b--');
end

for b = 1:size(betas,2)
    l2 = plot3(Ay_mat(b,:),Cn_mat(b,:),Ax_mat(b,:),'r');
end

legend([l1,l2],{'Constant Steer','Constant Vehicle Slip'})
hold off
  
