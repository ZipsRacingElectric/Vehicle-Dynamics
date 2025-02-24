%{
%% Overview:
This block computes the slip angle in radians for each tire.

y_dot: vehicle cg velocity along vehicle y axis (m/s)
x_dot: vehicle cg velocity along vehicle x axis (m/s)

r: vehicle yaw rate (rad/s)

delta_fl: front left wheel steering angle (rad)
delta_fr: front right wheel steering angle (rad)
delta_rl: rear left wheel steering angle (rad)
delta_rr: rear right wheel steering angle (rad)

a: distance of front wheels to cg viewed from the top
b: distance of rear wheels to cg viewed from the top

tf: front track width
tr: rear track width

alpha_fl: front wheel slip angle
alpha_fr: front wheel slip angle
alpha_rl: rear wheel slip angle
alpha_rr: rear wheel slip angle
%}

function [alpha_fl, alpha_fr, alpha_rl, alpha_rr] = slip_angle(y_dot, x_dot, r, delta_fl, delta_fr, delta_rl, delta_rr, a, b, tf, tr)

    tol = 0.01; % tolerance for edge case conditions

    % Calculate the vehicle speed magnitude
    speed = sqrt(x_dot^2 + y_dot^2);
    
    % If speed is nearly zero, define all slip angles as zero
    if speed < tol
        alpha_fl = 0;
        alpha_fr = 0;
        alpha_rl = 0;
        alpha_rr = 0;
        return;
    end

    % Front left tire
    denom_fl = x_dot + r * tf/2;
    if abs(denom_fl) < tol
        alpha_fl = 0;
    else
        alpha_fl = atan((y_dot + a * r) / denom_fl) - delta_fl;
    end

    % Front right tire
    denom_fr = x_dot - r * tf/2;
    if abs(denom_fr) < tol
        alpha_fr = 0;
    else
        alpha_fr = atan((y_dot + a * r) / denom_fr) - delta_fr;
    end

    % Rear left tire
    denom_rl = x_dot + r * tr/2;
    if abs(denom_rl) < tol
        alpha_rl = 0;
    else
        alpha_rl = atan((y_dot - b * r) / denom_rl) - delta_rl;
    end

    % Rear right tire
    denom_rr = x_dot - r * tr/2;
    if abs(denom_rr) < tol
        alpha_rr = 0;
    else
        alpha_rr = atan((y_dot - b * r) / denom_rr) - delta_rr;
    end
end
