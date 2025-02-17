%{
%% Overview:
This script generates CSAPS splines functions for combined cornering, as a
function of 4 independent variables: FZ, SL, SA, and IA. Pressure is held
constant.

When you run this script, a lookup table of the CSAPS spline is saved in
the same directory as d2704_7in_csaps.mat

Why CSAPS over a more traditional model? Pacejeka models dont like assymetric tire data
(big problem for D2704), don't work well with Mz, and can only get 2 of the
3 characteristics of the tire curve correctly (linear region, transition region, and saturation region).
A CSAPS function results in a much better fit. Additionally, CSAPS splines 
can be generated on any data that has been re-arranged. This means you can 
create inverse tire models easilly, where for example you input a desired long. force and recieve a
target slip ratio to drive a wheel at (very useful for torque vectoring).

It has been modified from Bill Cobb's matlab script (availible on the TTC forum) 
for use with round 9 data.

%% Notes:
- load the formatted_data.mat file, generated using the data_formatting.m
script in this same directory
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- SAE sign convention 
- formatted_data expected columns(in order): SL, SA, IA, P, Fz, Fy, Fx, Mz, Mx

%% TODO:
- Mz lookup table
- re-order the pure and combined datasets so they eventually can be
combined in to one large lookup table per tire outputs

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
clear fmdata;

%% load combined cornering fitted data points
load("../tire_data/d2704/d2704_7in_formatted_combined_data.mat");

fmdata = fmdata_combined;

%% select data that holds P = 8 psi constant (eliminate a dimension to simiplify lookup table)
srs = unique(round(fmdata(:,1), 2));
nsrs = length(srs); % number of distinct slip ratios
slips = unique(round(fmdata(:,2)))';
nslips = length(slips); % figure out number of distinct slip angles
incls = unique(round(fmdata(:,3)))';
nincls = length(incls); % figure out number of distinct inclination angles
press = unique(round(fmdata(:,4)))';
npress = length(press); % figure out number of distinct pressures

p8 = find(fmdata(:,4) == 8); % find array indicies for 8 psi
fmdata_p8 = fmdata(p8,:); % selects all data corresponding to 8 psi

%% Next, we transpose the arrays to make our spline functions happy:
% first take the sorted Fz column data and reshape it into individual columns. Each column
% represents a distint slip angle, slip ratio, and inclination angle combo.
% This creates an array where all the measured normal loads are along rows,
% with each row for a distinct normal load.
loads = reshape(fmdata_p8(:,5),[],(nslips * nincls * nsrs));

% take the average of each row to come up with unique loads. Also
% transpose into a row vector.
loads = mean(loads,2)';
nloads = length(loads);

% reshape the Fy data into 5D load, slip ratio, slip angle, and IA dimensions representive of the
% surface csaps will try and fit
% A 3D surface representation might look like: Fy vs Fz, SL, for
% a fixed SL, IA, and P
fy_p8 = reshape(fmdata_p8(:,6), nloads, nslips, nincls, nsrs);


% permutate (switch) the rows and columns such that we get:
% rows: SA sweep
% colums: Fz sweep
% vertical dimension: IA
% 4th dimension: SL sweep
fy_p8 = permute(fy_p8, [2 1 3 4]);

%% FY Surface Fit (P = 8 psi)
Fy_from_sa_fz_ia_sl = csaps({slips, loads, incls, srs}, fy_p8);

sl_value = 0.2;
ia_value = 0;
Fy_surface = fnval(Fy_from_sa_fz_ia_sl, {slips, loads, ia_value, sl_value});

figure
surf(loads, slips, Fy_surface)
xlabel('Normal Load (N)');
ylabel('Slip Angle (deg)');
zlabel('Lateral Force (N)')
title(sprintf('CSAPS Fit -  Lateral Force vs. Slip Angle, Load (SL = %d, IA = %d)', sl_value, ia_value));
view(45,45)
colorbar

%% Generate higher resolution lookup table data
% evaluate the csaps function at the resolution we want in each dimension
% store the data in a structure we can import into a Simulink n-D lookup block
% these ranges should match the pure cornering lookup table so they can be
% combined later

sl_min = -0.3;    % Min value for 'sl'
sl_max = 0.3;    % Max value for 'sl'
sl_res = 0.06;    % Resolution for 'sl'

sa_min = -13;    % Min value for 'sa'
sa_max = 13;     % Max value for 'sa'
sa_res = 0.25;    % Resolution for 'sa'

fz_min = -1600;  % Min value for 'fz'
fz_max = 0;   % Max value for 'fz'
fz_res = 25;    % Resolution for 'fz'

ia_min = 0;      % Min value for 'ia'
ia_max = 4;      % Max value for 'ia'
ia_res = 0.25;    % Resolution for 'ia'

% create grid of all possible input combinations
sl = sl_min:sl_res:sl_max;
sa = sa_min:sa_res:sa_max;
fz = fz_min:fz_res:fz_max;
ia = ia_min:ia_res:ia_max;

[SA, FZ, IA, SL] = ndgrid(sa, fz, ia, sl);

fy = zeros(size(SL));  % Preallocate the output matrix, size is the total number of possible input combinations

% Loop through each combination of inputs and evaluate the function
for i = 1:numel(SL)
    fy(i) = fnval(Fy_from_sa_fz_ia_sl, {SA(i), FZ(i), IA(i), SL(i)});
end

% validate lookup table data before saving
% plots the same surface we evaluated the csaps spline for above
[~, ia_index] = min(abs(ia - ia_value));  % Find the index of the closest 'ia_value'
[~, sl_index] = min(abs(sl - sl_value));  % Find the index of the closest 'ia_value'
fy_fixed_ia = fy(:, :, ia_index, sl_index)';  % Extract 2D matrix for fixed 'ia and sl, and transpose to match dimensions
[SA_grid, FZ_grid] = meshgrid(sa, fz);

% Using surf for a 3D surface plot
figure;
surf(FZ_grid, SA_grid, fy_fixed_ia);
xlabel('Normal Load (N)');
ylabel('Slip Angle (deg)');
zlabel('Lateral Force (N)');
title(sprintf('Lookup Table -  Lateral Force vs. Slip Angle, Load (SL = %d, IA = %d)', sl_value, ia_value));
view(45,45)
colorbar

% save both spline and lookup data in .mat file for use in simulink
save("../tire_data/d2704/d2704_7in_csaps_combined.mat", 'Fy_from_sa_fz_ia_sl', 'sl', 'sa', 'fz', 'ia', 'fy');