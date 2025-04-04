%{
%% Overview
This script uses the symbolic toolbox to derive the QP problem matricies given
the objective and constraint definitions.

%% Details:

%% TODO:
- Vehicle dynamics for Ax and SL objective functions
- Define constraints
- generate C function for using QP matricies
%}


clc; clear all;

%% Vehicle Dynamics
syms Fy_fl Fy_fr Fy_rl Fy_rr Fx_fl Fx_fr Fx_rl Fx_rr real
syms delta_fl delta_fr delta_rl delta_rr real
syms lf lr tf tr R real
syms mass real

% Define tire distance vectors in vehicle coordinate system
d_fl = [lf; tf/2; 0]; % x, y, z
d_fr = [lf; -tf/2; 0];
d_rl = [-lr; tr/2; 0];
d_rr = [-lr; -tr/2; 0];

% Define tire force vectors in tire coordinate system
Ft_fl = [Fx_fl; Fy_fl; 0]; % x, y, z
Ft_fr = [Fx_fr; Fy_fr; 0];
Ft_rl = [Fx_rl; Fy_rl; 0];
Ft_rr = [Fx_rr; Fy_rr; 0];

% Solve for tire force vectors in the vehicle coordinate system
F_fl = tire_to_vehicle(Ft_fl, delta_fl); % x, y, z
F_fr = tire_to_vehicle(Ft_fr, delta_fr);
F_rl = tire_to_vehicle(Ft_rl, delta_rr);
F_rr = tire_to_vehicle(Ft_rr, delta_rr);

% Solve for individual contributions to Mz
Mz_fl = cross(d_fl, F_fl);
Mz_fr = cross(d_fr, F_fr);
Mz_rl = cross(d_rl, F_rl);
Mz_rr = cross(d_rr, F_rr);

% Solve for longitudinal forces of vehicle
Fx = F_fl(1) + F_fr(1) + F_rl(1) + F_rr(1);

%% Decision Variables
% We will define the matrix T later in the QP form section
syms T_fl T_fr T_rl T_rr real

%% Linearizing tire forces
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

% Solve for slip ratios for the 3rd objective function
s_fl_l = solve(eq_fl, s_fl);
s_fr_l = solve(eq_fr, s_fr);
s_rl_l = solve(eq_rl, s_rl);
s_rr_l = solve(eq_rr, s_rr);

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

Fx = subs(Fx, [Fx_fl, Fy_fl], [T_fl/R, Fy_fl_l]);
Fx = subs(Fx, [Fx_fr, Fy_fr], [T_fr/R, Fy_fr_l]);
Fx = subs(Fx, [Fx_rl, Fy_rl], [T_rl/R, Fy_rl_l]);
Fx = subs(Fx, [Fx_rr, Fy_rr], [T_rr/R, Fy_rr_l]);

% Note that this does not include tire Mz forces, this is done to simplify
% torque optimization, and they do not contribute significantly to vehicle
% Mz
Mz = Mz_fl(3) + Mz_fr(3) + Mz_rl(3) + Mz_rr(3);
ax = Fx / mass;
T = [T_fl;
     T_fr; 
     T_rl;
     T_rr];

%% Yaw Moment Objective Function
% our objective is f1 = (Mz_tv - Mz)^2 = (A1T - b1)^2
syms w1 alpha_1 Mz_tv real
A1 = sym('A1_', [1 4]);
b1 = sym('b', [1 1]);

% Convert Mz into matrix form Mz = A_mz*T + b_mz, then solve for A1 and b1
[A_mz, b_mz] = equationsToMatrix(Mz, T);
A1 = -A_mz;
b1 = b_mz - Mz_tv;

%% Longitudinal Objective Function
% f2 = (Ax_ref - Ax)^2 = (A2T - b2)^2
syms w2 alpha_2 ax_ref real
A2 = sym('A2_', [1 4]);
b2 = sym('b', [1 1]);

% Convert Ax into matrix form Ax = A_ax*T + b_ax, then solve for A2 and b2
[A_ax, b_ax] = equationsToMatrix(ax, T);
A2 = -A_ax;
b2 = b_ax - ax_ref;

%% Slip Ratio Objective Function
syms w3 alpha_3 real
% f3 = (A3T - b3)^2, so we convert it to matrix form.
% Note: f3 is already non-linear, so we cannot directly use
% equationsToMatrix()
% Instead we solve for the coefficients of T_fl, T_fr, etc.... and then build
% the matrix form manually

% Use coeffs() to extract the coefficient of T_i and the constant term for
% each slip angle
[cf_fl, terms_fl] = coeffs(s_fl_l, T_fl);
m_fl = sym(0); b_fl = sym(0);
for k = 1:length(terms_fl)
    if isequal(terms_fl(k), T_fl)
        m_fl = cf_fl(k);
    elseif isequal(terms_fl(k), sym(1))
        b_fl = cf_fl(k);
    end
end

[cf_fr, terms_fr] = coeffs(s_fr_l, T_fr);
m_fr = sym(0); b_fr = sym(0);
for k = 1:length(terms_fr)
    if isequal(terms_fr(k), T_fr)
        m_fr = cf_fr(k);
    elseif isequal(terms_fr(k), sym(1))
        b_fr = cf_fr(k);
    end
end

[cf_rl, terms_rl] = coeffs(s_rl_l, T_rl);
m_rl = sym(0); b_rl = sym(0);
for k = 1:length(terms_rl)
    if isequal(terms_rl(k), T_rl)
        m_rl = cf_rl(k);
    elseif isequal(terms_rl(k), sym(1))
        b_rl = cf_rl(k);
    end
end

[cf_rr, terms_rr] = coeffs(s_rr_l, T_rr);
m_rr = sym(0); b_rr = sym(0);
for k = 1:length(terms_rr)
    if isequal(terms_rr(k), T_rr)
        m_rr = cf_rr(k);
    elseif isequal(terms_rr(k), sym(1))
        b_rr = cf_rr(k);
    end
end

% Form the A3 and b3 in f3 = (A3*T - b3)^2
% We want A3*T - b = [m_fl*T_fl + b_fl; m_fr*T_fr + b_fr; m_rl*T_rl + b_rl; m_rr*T_rr + b_rr]
A3 = diag([m_fl, m_fr, m_rl, m_rr]);
b3 = -[b_fl; b_fr; b_rl; b_rr];

% Now form f3 in matrix form and compare with the sum of squares form.
f3_matrixForm = (A3*T - b3).'*(A3*T - b3);
f3_original = s_fl_l^2 + s_fr_l^2 + s_rl_l^2 + s_rr_l^2;

% Verify the equivalence:
difference = simplify(f3_original - f3_matrixForm);
disp('Difference (should be 0):');
pretty(difference)

%% Constraints

%% QP Problem Form
% QP form is 1/2 * T^T * P * T + q^T * T. Also, P = 2*Q. We derive this manually
% because matlab sucks

% Compute Q and P: the quadratic coefficient matrix
Q = (w1/alpha_1)*(A1.' * A1) + (w2/alpha_2)*(A2.' * A2) + (w3/alpha_3)*(A3.' * A3);
P = simplify(2 * Q)

% Verify P is symmetric (aka everything worked as expected)
diffP = simplify(P - P.');
disp('P should be symmetric (output should be 0):');
disp(diffP)

% Compute q: the linear coefficient vector (as a column vector)
q = simplify(-2 * ((w1/alpha_1)*(A1.' * b1) + (w2/alpha_2)*(A2.' * b2) + (w3/alpha_3)*(A3.' * b3)))

%% Save QP Matricies
save('QP_matrices.mat', 'P', 'q');

%% Helper Functions
function force_v = tire_to_vehicle(force_t, delta)
    rotation = [cos(delta), -sin(delta),  0;
                sin(delta),  cos(delta),  0;
                         0,          0,   1];
    force_v = rotation * force_t;
end