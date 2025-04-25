%{

%% Overview
This script graphs lateral acceleration as a function of the understeer
gradient for fixed steering wheel angles. It is intended to show that the ideal
understeer gradient is Ku = 0 for any vehicle.


%% Details:
- Linear single-track dynamic model is how the understeer gradient is
derived from.

%% TODO:

%}

clear;  clc;  close all

% ensure vehicle object in in matlab path
directory = fileparts(which(mfilename)); % Get the directory of the current script
parentDirectory = fileparts(directory); % Get the parent directory of the current script
addpath(genpath(parentDirectory));      % Add the parent directory and all its subdirectories to the MATLAB path

% load vehicle data
zr25 = vehicle("../vehicle_data/zr25_data.xlsx");

L        = zr25.wheelbase;              % wheel-base  [m]
R        = 8;                          % path radius [m]  (fixed)
steerDeg = [5 10 15 18.7];                % road-wheel steer angles to plot [deg]
mu       = 1.5;                         % tyre-road friction coefficient  (≈1.05 on good asphalt)
g        = zr25.g;                      % gravity [m/s^2]
KuRange  = linspace(-0.02,0.10,4001);   % under-/over-steer gradient sweep  [rad/(m/s^2)]
KuEps    = 1e-5;                        % exclude a tiny region around Ku = 0

deltaKin = L/R;                         % pure‐kinematic steer [rad]
deltaDyn = deg2rad(steerDeg(:));        % column vector of dynamic steers [rad]

ay = nan(numel(deltaDyn),numel(KuRange));

for k = 1:numel(deltaDyn)
    C = deltaDyn(k) - deltaKin;         % numerator from derivation
    KuValid = abs(KuRange) > KuEps;     % avoid division by ~0
    ay_k    = C ./ KuRange;             % eqn  ay = (deltaDyn-L/R)/Ku
    ay_k(~KuValid) = NaN;               % punch a hole at Ku ≈ 0
    ay_k(ay_k < 0) = NaN;               % ignore negative (physically unusable) values
    ay_k(ay_k > mu*g) = mu*g;           % impose grip limit (cap)
    ay(k,:) = ay_k;
end

% Plot ay vs Ku with a maximum grip limit
figure;  hold on;  box on
colors = lines(numel(deltaDyn));
for k = 1:numel(deltaDyn)
    plot(KuRange, ay(k,:)/g, 'LineWidth',1.5, 'Color',colors(k,:), ...
         'DisplayName', sprintf('\\delta = %.0f°', steerDeg(k)));
end
yline(mu,'k--','LineWidth',1.2,'DisplayName','tire grip limit');

xlabel('Understeer gradient  K_u  [rad/(m/s^2)]');
ylabel('Lateral acceleration  A_y [g]');
title(sprintf('a_y vs K_u for fixed radius = %.0f m',R));
legend('Location','northeast');
grid on;  xlim([min(KuRange) max(KuRange)]);
ylim([0 1.2*mu]);