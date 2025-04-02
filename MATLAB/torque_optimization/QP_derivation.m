%{

This script uses the symbolic toolbox to derive the QP problem matricies given
the objective and constraint definitions.

%}

clc; clear all;

%% Vehicle Dynamics
syms Fy_fl Fy_fr Fy_rl Fy_rr Fx_fl Fx_fr Fx_rl Fx_rr real
syms delta_fl delta_fr delta_rl delta_rr real
syms lf lr tf tr R real

% Define tire distance vectors in vehicle coordinate system
d_fl = [lf; tf/2; 0];
d_fr = [lf; -tf/2; 0];
d_rl = [-lr; tr/2; 0];
d_rr = [-lr; -tr/2; 0];

% Define tire force vectors in tire coordinate system
Ft_fl = [Fx_fl; Fy_fl; 0];
Ft_fr = [Fx_fr; Fy_fr; 0];
Ft_rl = [Fx_rl; Fy_rl; 0];
Ft_rr = [Fx_rr; Fy_rr; 0];

% Solve for tire force vectors in the vehicle coordinate system
F_fl = tire_to_vehicle(Ft_fl, delta_fl);
F_fr = tire_to_vehicle(Ft_fr, delta_fr);
F_rl = tire_to_vehicle(Ft_rl, delta_rr);
F_rr = tire_to_vehicle(Ft_rr, delta_rr);

% Solve for individual contributions to Mz
Mz_fl = cross(d_fl, F_fl);
Mz_fr = cross(d_fr, F_fr);
Mz_rl = cross(d_rl, F_rl);
Mz_rr = cross(d_rr, F_rr);

%% Decision Variables
% We will define the matrix T later in the QP form section
syms T_fl T_fr T_rl T_rr real

%% Linearizing Mz contributions
% When linearizing, we define Fy and Fx as linear functions around their
% operating point:
syms Fx_fl_0 Fx_fr_0 Fx_rl_0 Fx_rr_0 Fy_fl_0 Fy_fr_0 Fy_rl_0 Fy_rr_0 real
syms s_fl_0 s_fr_0 s_rl_0 s_rr_0 s_fl s_fr s_rl s_rr real
syms dFx_fl_ds dFx_fr_ds dFx_rl_ds dFx_rr_ds dFy_fl_ds dFy_fr_ds dFy_rl_ds dFy_rr_ds real

% Define linear versions in Point-slope form:
Fx_fl_l = dFx_fl_ds * (s_fl - s_fl_0) + Fx_fl_0;
Fx_fr_l = dFx_fr_ds * (s_fr - s_fr_0) + Fx_fr_0;
Fx_rl_l = dFx_rl_ds * (s_rl - s_rl_0) + Fx_rl_0;
Fx_rr_l = dFx_rr_ds * (s_rr - s_rr_0) + Fx_rr_0;

Fy_fl_l = dFy_fl_ds * (s_fl - s_fl_0) + Fy_fl_0;
Fy_fr_l = dFy_fr_ds * (s_fr - s_fr_0) + Fy_fr_0;
Fy_rl_l = dFy_rl_ds * (s_rl - s_rl_0) + Fy_rl_0;
Fy_rr_l = dFy_rr_ds * (s_rr - s_rr_0) + Fy_rr_0;

% Solve for Fy as functions of T
eq_fl = Fx_fl_l == T_fl / R;
eq_fr = Fx_fr_l == T_fr / R;
eq_rl = Fx_rl_l == T_rl / R;
eq_rr = Fx_rr_l == T_rr / R;

sol_fl = simplify(solve(eq_fl, s_fl) - s_fl_0);
sol_fr = simplify(solve(eq_fr, s_fr) - s_fr_0);
sol_rl = simplify(solve(eq_rl, s_rl) - s_rl_0);
sol_rr = simplify(solve(eq_rr, s_rr) - s_rr_0);

Fy_fl_l = subs(Fy_fl_l, (s_fl - s_fl_0), sol_fl);
Fy_fr_l = subs(Fy_fr_l, (s_fr - s_fr_0), sol_fr);
Fy_rl_l = subs(Fy_rl_l, (s_rl - s_rl_0), sol_rl);
Fy_rr_l = subs(Fy_rr_l, (s_rr - s_rr_0), sol_rr);

% Sub in our linear functions for Fx and Fy
% Note Fx can directly be substituted with Fx_i = Ti / R
Mz_fl = subs(Mz_fl, [Fx_fl, Fy_fl], [T_fl/R, Fy_fl_l]);
Mz_fr = subs(Mz_fr, [Fx_fr, Fy_fr], [T_fr/R, Fy_fr_l]);
Mz_rl = subs(Mz_rl, [Fx_rl, Fy_rl], [T_rl/R, Fy_rl_l]);
Mz_rr = subs(Mz_rr, [Fx_rr, Fy_rr], [T_rr/R, Fy_rr_l]);

% Note that this does not include tire Mz forces, this is done to simplify
% torque optimization, and they do not contribute significantly to vehicle
% Mz
Mz = Mz_fl(3) + Mz_fr(3) + Mz_rl(3) + Mz_rr(3);

%% Objective Functions
T = sym('T', [4 1]);
% Note: all have the form f = (AT - b)^2, which ensures that the P matrix
% is always positive semi-definite

% Yaw Moment Objective Function
% our objective is f1 = (Mz_tv - Mz)^2 = (A1T - b1)^2
syms w1 alpha_1 Mz_tv real
A1 = sym('A1_', [1 4]);
b1 = sym('b', [1 1]);

% Matrix form for Mz:
[A_mz, b_mz] = equationsToMatrix(Mz, [T_fl, T_fr, T_rl, T_rr]);
A1 = -A_mz;
b1 = b_mz - Mz_tv;

% We expand this manually because matlab sucks. Note, we left off the
% constant portion because it does not change the optimization at all
f1 = T.' * (A1.' * A1) * T - 2 * T.' * (A1.' * b1);

% Longitudinal Objective Function
% f2 = (Ax_ref - Ax)^2 = (A2T - b2)^2
syms w2 alpha_2 b2 real
A2 = sym('A2_', [1 4]);

% TODO
%A2 = 
%b2 = 

f2 = T.' * (A2.' * A2) * T - 2 * T.' * (A2.' * b2);

% Longitudinal Objective Function
% f2 = s_fl^2 + s_fr^2 + s_rl^2 + s_rr^2 = (A3T - b3)^2
syms w3 alpha_3 b3 real
A3 = sym('A3_', [1 4]);

% TODO
%A3 = 
%b3 = 

f3 = T.' * (A3.' * A3) * T - 2 * T.' * (A3.' * b3);

%% Constraints

%% QP Problem Form
% QP form is 1/2 * T^T * P * T + q^T * T. Also, P = 2*Q. We also derive this manually
% because matlab sucks

% Compute Q: the quadratic coefficient matrix
Q = (w1/alpha_1)*(A1.' * A1) + (w2/alpha_2)*(A2.' * A2) + (w3/alpha_3)*(A3.' * A3)

% Compute q: the linear coefficient vector (as a column vector)
q = -2 * ( (w1/alpha_1)*(A1.' * b1) + (w2/alpha_2)*(A2.' * b2) + (w3/alpha_3)*(A3.' * b3) );

% Now, populate the Q and q matricies the cost function coefficients
% TODO: VD for f2 and f3 so we can sub those in
%Q = subs(A1(1), -A_mz(1))
%q = subs(b1, (b_mz - Mz_tv))

%% Helper Functions
function force_v = tire_to_vehicle(force_t, delta)
    rotation = [cos(delta), -sin(delta),  0;
                sin(delta),  cos(delta),  0;
                         0,          0,   1];
    force_v = rotation * force_t;
end
