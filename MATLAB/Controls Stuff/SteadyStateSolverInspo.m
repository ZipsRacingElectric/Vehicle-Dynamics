%% 4WD EV Steady-State Solver (Vehicle-Level Outputs)
% Inputs:
%   Vehicle - parameter object/struct
%   SAS   - front steering angle [rad]
%   Vx      - longitudinal speed [m/s]


Vehicle = zr25
SteerAng = 2
Vx = 20 %m/s
[Vy_steady, r_steady, ay_steady] = steadyStateVehicleLevel(Vehicle, SteerAng, Vx)



%% Functions

function [Vy_steady, r_steady, ay_steady] = steadyStateVehicleLevel(Vehicle, SteerAng, Vx)

    % Extract vehicle parameters
    m   = Vehicle.mass;
    Izz = Vehicle.Izz;
    lf  = Vehicle.lf;
    lr  = Vehicle.lr;
    Cf  = Vehicle.Cf; % total front cornering stiffness
    Cr  = Vehicle.Cr; % total rear cornering stiffness

    % Initial guess [Vy, r]
    x0 = [0; 0];

    % Residual function for fsolve
    fun = @(x) lateralYawResidual(x, Vx, SteerAng, m, Izz, lf, lr, Cf, Cr);

    options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-8);
    sol = fsolve(fun, x0, options);

    % Extract steady-state values
    Vy_steady = sol(1);
    r_steady  = sol(2);

    % Compute lateral acceleration at CG
    ay_steady = Vx * r_steady;  % planar model approximation

end

%% --- Residual function for lateral & yaw balance ---
function res = lateralYawResidual(x, Vx, delta, m, Izz, lf, lr, Cf, Cr)
    Vy = x(1);
    r  = x(2);

    % Slip angles (front/rear) for linear bicycle approximation
    alpha_f = delta - (Vy + lf * r)/Vx;
    alpha_r = - (Vy - lr * r)/Vx;

    % Lateral forces (linear tires)
    Fyf = Cf * alpha_f;
    Fyr = Cr * alpha_r;

    % Residuals: sum of lateral forces and yaw moments
    res1 = Fyf + Fyr;                  % lateral force balance (m*Vy_dot ~ 0)
    res2 = lf*Fyf - lr*Fyr;            % yaw moment balance (Izz*r_dot ~ 0)

    res = [res1; res2];
end
