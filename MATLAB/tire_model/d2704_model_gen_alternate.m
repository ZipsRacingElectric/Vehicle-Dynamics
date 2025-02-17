%{
%% Overview:
This script attemps to solve some issues combining the two datasets by
using radius basis function interpolation instead of a cubic
spline interpolation. This has the advantage of interpolating on scattered
data outside the convex radius.

%% Notes:
- This script inputs the formatted datasets created with the formatting
scripts.
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- SAE sign convention 
- formatted_data expected columns(in order): SL, SA, IA, P, Fz, Fy, Fx, Mz, Mx
- The output matricies are in the following dimension order: SL, SA, Fz

%% TODO:
- sweep a local region of sigma value and graph the resulting error
- The Fx fit when SL = 0 is not great for some normal loads. Need to figure
out how to populate this with data points as we know it is = to zero,
however Fy, Mz are not zero, so this isn't straightforward

%% Data and units:
FX: longituidinal force, N
FY: lateral force, N
FZ: normal load, N
IA: inclination angle, degrees
MX: overturning moment, N-m
MZ: aligning torque, N-m
P: tire pressure, psi
SA: slip angle, degrees
SL: slip ratio, unitless
%}

clear; close all; clc;

%% Add RBF scripts to path
addpath('rbfinterp_v1.2'); % Change to your actual folder path


%% load fitted data points
% Format: SL, SA, IA, P, Fz, Fy, Fx, Mz, Mx
load("../tire_data/d2704/d2704_7in_formatted_combined_data.mat");
load("../tire_data/d2704/d2704_7in_formatted_data.mat");

%% Combine data and simplify
% First delete SL = 0 from fmdata_combined since we want to use the pure
% conering data for this condition.
indexes = fmdata_combined(:,1) ~= 0;
fmdata_combined = fmdata_combined(indexes,:);

% Combined datasets
fmdata = [fmdata_pure; fmdata_combined];

% Select P = 8 psi and IA = 2 deg
fmdata = fmdata(find(fmdata(:,4) == 8 & fmdata(:,3) == 2),:);

% Remove IA and P columns since we wont need them anymore. Remove Mx
% column since we wont be using that data and want the simplelest thing
% possible
fmdata(:, [3,4, 9]) = [];

% Sort by Slip Ratio, Slip Angle, and finally Fz
% fmdata is now [SL, SA, Fz, Fy, Fx, Mz]
fmdata = sortrows(fmdata, [1,2,3]);

%% Normalization of input data
% Extract independent variables (inputs)
SL = fmdata(:,1);  % Slip Ratio
SA = fmdata(:,2);  % Slip Angle
Fz = fmdata(:,3);  % Normal Load

% Normalize
SL_norm = (SL - min(SL)) ./ (max(SL) - min(SL));
SA_norm = (SA - min(SA)) ./ (max(SA) - min(SA));
Fz_norm = (Fz - min(Fz)) ./ (max(Fz) - min(Fz));

%% Radial Base Function Interpolation
% Extract dependent variables (outputs)
Fy = fmdata(:,4);  % Lateral Force
Fx = fmdata(:,5);  % Longitudinal Force
Mz = fmdata(:,6);  % Aligning Moment

% Combine independent variables into a single matrix (of sample points)
X = [SL_norm, SA_norm, Fz_norm]';

% adjust these to adjust model lfit
sigma = 0.03;
smooth_param = 0.7;

% Train RBF models for Fy, Fx, and Mz
rbf_Fy = rbfcreate(X, Fy', 'RBFFunction', 'multiquadric', 'RBFConstant', sigma, 'RBFSmooth', smooth_param);
rbf_Fx = rbfcreate(X, Fx', 'RBFFunction', 'multiquadric', 'RBFConstant', sigma, 'RBFSmooth', smooth_param);
rbf_Mz = rbfcreate(X, Mz', 'RBFFunction', 'multiquadric', 'RBFConstant', sigma, 'RBFSmooth', smooth_param);

% Define query point for interpolation/extrapolation (must be normalized
SL_q = (0.05 - min(SL)) ./ (max(SL) - min(SL));  % Example Slip Ratio
SA_q = (5 - min(SA)) ./ (max(SA) - min(SA));     % Example Slip Angle
Fz_q = (2000 - min(Fz)) ./ (max(Fz) - min(Fz));  % Example Normal Load
query_point = [SL_q; SA_q; Fz_q];

% Evaluate RBF models at the query point
Fy_q = rbfinterp(query_point, rbf_Fy);
Fx_q = rbfinterp(query_point, rbf_Fx);
Mz_q = rbfinterp(query_point, rbf_Mz);

% Display results
fprintf('RBF at sample query point:\n Fy = %.2f, Fx = %.2f, Mz = %.2f\n', Fy_q, Fx_q, Mz_q);

%% Generate Lookup Table Data
sl_min = -0.3;    % Min value for 'sl'
sl_max = 0.3;    % Max value for 'sl'
sl_res = 0.02;    % Resolution for 'sl'

sa_min = -13;    % Min value for 'sa'
sa_max = 13;     % Max value for 'sa'
sa_res = 0.25;    % Resolution for 'sa'

fz_min = -1600;  % Min value for 'fz'
fz_max = 0;   % Max value for 'fz'
fz_res = 25;    % Resolution for 'fz'

% Define lookup table grid ranges
sl_vals = sl_min:sl_res:sl_max;
sa_vals = sa_min:sa_res:sa_max;
fz_vals = fz_min:fz_res:fz_max;

% Create a meshgrid of SL, SA, and Fz
[SL_grid, SA_grid, Fz_grid] = ndgrid(sl_vals, sa_vals, fz_vals);

% Normalize grid inputs (must match RBF normalization)
SL_grid_norm = (SL_grid - min(SL)) / (max(SL) - min(SL));
SA_grid_norm = (SA_grid - min(SA)) / (max(SA) - min(SA));
Fz_grid_norm = (Fz_grid - min(Fz)) / (max(Fz) - min(Fz));

% Reshape grid to pass as query points
query_points_norm = [SL_grid_norm(:), SA_grid_norm(:), Fz_grid_norm(:)]';

% Evaluate RBF models for each grid point
Fy_grid = rbfinterp(query_points_norm, rbf_Fy);
Fx_grid = rbfinterp(query_points_norm, rbf_Fx);
Mz_grid = rbfinterp(query_points_norm, rbf_Mz);

% Reshape interpolated values to match grid size
Fy_grid = reshape(Fy_grid, size(SL_grid));
Fx_grid = reshape(Fx_grid, size(SL_grid));
Mz_grid = reshape(Mz_grid, size(SL_grid));

%% Fy vs SA, Fz, SL = 0
sl_value = 0;
[~, sl_index] = min(abs(sl_vals - sl_value));  % Find the index of the closest 'ia_value'
Fy_fixed_sl = squeeze(Fy_grid(sl_index, :, :));

% Using surf for a 3D surface plot
figure;
surf(fz_vals, sa_vals, Fy_fixed_sl);
hold on;
xlabel('Normal Load (N)');
ylabel('Slip Angle (deg)');
zlabel('Lateral Force (N)');
title(sprintf('Lookup Table -  Lateral Force vs. Slip Angle, Load (SL = %d, IA = %d), sigma = %d', sl_value, 2, sigma));
view(45,45)
colorbar

% Plot the original sample points
sl_idx = SL == sl_value; % Find samples close to sl_index
scatter3(Fz(sl_idx), SA(sl_idx), Fy(sl_idx), ...
    50, 'r', 'filled', 'MarkerEdgeColor', 'k'); % Red markers with black edges
grid on;
legend('Interpolated Surface', 'Original Sample Points');
hold off;

%% Save the data
save("../tire_data/d2704/d2704_7in_rbf_full_model.mat", "sl_vals", "sa_vals", "fz_vals", "Fy_grid", "Fx_grid", "Mz_grid");

%% Plot a bunch of other graphs

% Define SL, SA, and Fz values for fixed conditions
sl_values = [0, 0.1, 0.2, 0.3];      % Slip Ratios for Fy vs SA, Fz
sa_values = [0, 6, 12];              % Slip Angles for Fy vs SL, Fz
fz_values = [-300, -600, -1200];     % Normal Loads for Fy vs SL, SA

%% **1. Fy vs SA, Fz for SL = [0, 0.1, 0.2, 0.3]**
for sl_value = sl_values
    [~, sl_index] = min(abs(sl_vals - sl_value)); % Find closest index
    Fy_fixed_sl = squeeze(Fy_grid(sl_index, :, :)); % Extract data slice
    
    % Find original sample points for SL ≈ sl_value
    sl_idx = abs(SL - sl_value) < 1e-3;
    
    % Plot surface
    plot_surface(fz_vals, sa_vals, Fy_fixed_sl, Fz(sl_idx), SA(sl_idx), Fy(sl_idx), ...
        'Normal Load (N)', 'Slip Angle (deg)', 'Lateral Force Fy (N)', ...
        sprintf('Fy vs SA, Fz (SL = %.1f)', sl_value));
end

%% **2. Fx vs SA, Fz for SL = [0, 0.1, 0.2, 0.3]**
for sl_value = sl_values
    [~, sl_index] = min(abs(sl_vals - sl_value));
    Fx_fixed_sl = squeeze(Fx_grid(sl_index, :, :));
    sl_idx = abs(SL - sl_value) < 1e-3;
    
    plot_surface(fz_vals, sa_vals, Fx_fixed_sl, Fz(sl_idx), SA(sl_idx), Fx(sl_idx), ...
        'Normal Load (N)', 'Slip Angle (deg)', 'Longitudinal Force Fx (N)', ...
        sprintf('Fx vs SA, Fz (SL = %.1f)', sl_value));
end

%% **3. Mz vs SA, Fz for SL = [0, 0.1, 0.2, 0.3]**
for sl_value = sl_values
    [~, sl_index] = min(abs(sl_vals - sl_value));
    Mz_fixed_sl = squeeze(Mz_grid(sl_index, :, :));
    sl_idx = abs(SL - sl_value) < 1e-3;
    
    plot_surface(fz_vals, sa_vals, Mz_fixed_sl, Fz(sl_idx), SA(sl_idx), Mz(sl_idx), ...
        'Normal Load (N)', 'Slip Angle (deg)', 'Aligning Moment Mz (N·m)', ...
        sprintf('Mz vs SA, Fz (SL = %.1f)', sl_value));
end

%% **4. Fy vs SL, Fz for SA = [0, 6, 12]**
for sa_value = sa_values
    [~, sa_index] = min(abs(sa_vals - sa_value));
    Fy_fixed_sa = squeeze(Fy_grid(:, sa_index, :));
    sa_idx = abs(SA - sa_value) < 1e-3;
    
    plot_surface(fz_vals, sl_vals, Fy_fixed_sa, Fz(sa_idx), SL(sa_idx), Fy(sa_idx), ...
        'Normal Load (N)', 'Slip Ratio', 'Lateral Force Fy (N)', ...
        sprintf('Fy vs SL, Fz (SA = %d)', sa_value));
end

%% **5. Fx vs SL, Fz for SA = [0, 6, 12]**
for sa_value = sa_values
    [~, sa_index] = min(abs(sa_vals - sa_value));
    Fx_fixed_sa = squeeze(Fx_grid(:, sa_index, :));
    sa_idx = abs(SA - sa_value) < 1e-3;
    
    plot_surface(fz_vals, sl_vals, Fx_fixed_sa, Fz(sa_idx), SL(sa_idx), Fx(sa_idx), ...
        'Normal Load (N)', 'Slip Ratio', 'Longitudinal Force Fx (N)', ...
        sprintf('Fx vs SL, Fz (SA = %d)', sa_value));
end

%% **6. Mz vs SL, Fz for SA = [0, 6, 12]**
for sa_value = sa_values
    [~, sa_index] = min(abs(sa_vals - sa_value));
    Mz_fixed_sa = squeeze(Mz_grid(:, sa_index, :));
    sa_idx = abs(SA - sa_value) < 1e-3;
    
    plot_surface(fz_vals, sl_vals, Mz_fixed_sa, Fz(sa_idx), SL(sa_idx), Mz(sa_idx), ...
        'Normal Load (N)', 'Slip Ratio', 'Aligning Moment Mz (N·m)', ...
        sprintf('Mz vs SL, Fz (SA = %d)', sa_value));
end

%% **7. Fy vs SL, SA for Fz = [-300, -600, -1200]**
for fz_value = fz_values
    [~, fz_index] = min(abs(fz_vals - fz_value));
    Fy_fixed_fz = squeeze(Fy_grid(:, :, fz_index));
    fz_idx = abs(Fz - fz_value) < 1e-3;
    
    plot_surface(sa_vals, sl_vals, Fy_fixed_fz, SA(fz_idx), SL(fz_idx), Fy(fz_idx), ...
        'Slip Angle (deg)', 'Slip Ratio', 'Lateral Force Fy (N)', ...
        sprintf('Fy vs SL, SA (Fz = %d N)', fz_value));
end

%% **8. Fx vs SL, SA for Fz = [-300, -600, -1200]**
for fz_value = fz_values
    [~, fz_index] = min(abs(fz_vals - fz_value));
    Fx_fixed_fz = squeeze(Fx_grid(:, :, fz_index));
    fz_idx = abs(Fz - fz_value) < 1e-3;
    
    plot_surface(sa_vals, sl_vals, Fx_fixed_fz, SA(fz_idx), SL(fz_idx), Fx(fz_idx), ...
        'Slip Angle (deg)', 'Slip Ratio', 'Longitudinal Force Fx (N)', ...
        sprintf('Fx vs SL, SA (Fz = %d N)', fz_value));
end

% Define function for plotting interpolated surfaces and sample points
function plot_surface(x_vals, y_vals, z_vals, orig_x, orig_y, orig_z, x_label, y_label, z_label, title_str)
    figure;
    surf(x_vals, y_vals, z_vals, 'EdgeColor', 'none');
    hold on;
    scatter3(orig_x, orig_y, orig_z, 50, 'r', 'filled', 'MarkerEdgeColor', 'k'); % Overlay sample points
    xlabel(x_label);
    ylabel(y_label);
    zlabel(z_label);
    title(title_str);
    view(45, 45);
    colorbar;
    grid on;
    hold off;
end

