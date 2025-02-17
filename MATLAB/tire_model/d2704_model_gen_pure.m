%{
%% Overview:
This script generates CSAPS splines functions for pure cornering, as a
function of 3 independent variables: FZ, SA, and IA. Pressure is held
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
- The graphs overlay the CSAPS spline and the lookup table. They look like
they deviate, however its probably because the CSAPS spline is being
plotted for a different IA value. Need to debug plots

%% Data and units:
FX: longituidinal force, N
FY: lateral force, N
FZ: normal load, N
IA: inclination angle, degrees
MX: overturning moment, N-m
MZ: aligning torque, N-m
P: tire pressure, psi
SA: slip angle, degrees
%}

clear; close all; clc;
clear fmdata;

%% load pure cornering fitted data points
load("../tire_data/d2704/d2704_7in_formatted_data.mat");

fmdata = fmdata_pure;

%% select data that holds P = 8 psi constant, can also assume SL = 0
slips = unique(round(fmdata(:,2)))';
nslips = length(slips); % figure out number of distinct slip angles
incls = unique(round(fmdata(:,3)))';
nincls = length(incls); % figure out number of distinct inclination angles
press = unique(round(fmdata(:,4)))';
npress = length(press); % figure out number of distinct pressures

p8 = find(fmdata(:,4) == 8); % find array indexes for 8 psi
fmdata_p8 = fmdata(p8,:); % selects all data corresponding to 8 psi

%% Next, we transpose the arrays to make our spline functions happy:
% take the sorted Fz column data and reshape it into columns. Each column
% represents a distint slip angle and inclination angle combo
loads = reshape(fmdata_p8(:,5),[],(nslips*nincls));

% take the average of each row to come up with distinct loads. Also
% transpose into a row vector.
loads = mean(loads,2)';
nloads = length(loads);

% reshape the Fy data into 4D load, slip angle, and IA dimensions representive of the
% surface csaps will try and fit
% A 3D surface representation might look like: Fy fpr SA vs Fz, for
% a fixed SL, IA, and P
fy_p8 = reshape(fmdata_p8(:,6), nloads, nslips, nincls);


% permutate (switch) the rows and columns such that we get:
% rows: SA sweep
% colums: Fz sweep
% vertical dimension: IA
fy_p8 = permute(fy_p8, [2 1 3]);


%% FY Surface Fit (P = 8 psi)
Fy_from_sa_fz_ia = csaps({slips,loads,incls}, fy_p8);

sl_value = 0;
ia_value = 0;
Fy_surface = fnval(Fy_from_sa_fz_ia, {slips, loads, ia_value});

figure
surf(loads, slips, Fy_surface)
xlabel('Normal Load (N)');
ylabel('Slip Angle (deg)');
zlabel('Lateral Force (N)')
title(sprintf('CSAPS Fit -  Lateral Force vs. Slip Angle, Load (SL = %d, IA = %d)', sl_value, ia_value));
view(45,45)
colorbar

% evaluating the CSAPS tire model
fnval(Fy_from_sa_fz_ia,{6.3, -1200, 1.4})

%% Generate higher resolution lookup table data
% evaluate the csaps function at the resolution we want in each dimension
% store the data in a structure we can import into a Simulink n-D lookup block
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
sa = sa_min:sa_res:sa_max;
fz = fz_min:fz_res:fz_max;
ia = ia_min:ia_res:ia_max;

[SA, FZ, IA] = ndgrid(sa, fz, ia);

fy_pure = zeros(size(SA));  % Preallocate the output matrix

% Loop through each combination of inputs and evaluate the function
for i = 1:numel(SA)
    fy_pure(i) = fnval(Fy_from_sa_fz_ia, {SA(i), FZ(i), IA(i)});
end

% validate lookup table data before saving
% plots the same surface we evaluated the csaps spline for above
[~, ia_index] = min(abs(ia - ia_value));  % Find the index of the closest 'ia_value'
fy_fixed_ia = fy_pure(:, :, ia_index)';  % Extract 2D matrix for fixed 'ia_value' and transpose to match dimensions
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
save("../tire_data/d2704/d2704_7in_csaps.mat", 'Fy_from_sa_fz_ia', 'sa', 'fz', 'ia', 'fy_pure');

%{

%% Mz Surface Fit (P = 8 psi)
% reshape the data into 3D load, slip, and IA dimensions representive of the
% surface csaps will try and fit
mz_p8 = reshape(fmdata_p8(:,5),nloads,nslips,nincls); % column 5 is Mz data

% permutate (switch) the rows and columns such that we get:
% rows: SA sweep
% colums: Fz sweep
% vertical dimension: IA
mz_p8 = permute(mz_p8, [2 1 3]);

% generate csaps
Mz_from_sa_fz_ia = csaps({slips,loads,incls},mz_p8);

% note: when plotting the 4D function, it will give you a 3D surface at IA = 0
figure('Name','Overturning Moment (Mz) vs. Slip Angle, Vertical Load, Inclination Angle')
fnplt(Mz_from_sa_fz_ia)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force (N)')
title('Name','Overturning Moment (Mz) vs. Slip Angle, Vertical Load, IA = 0 deg')
view(45,45)

% evaluating the CSAPS tire model
% keep in mind the SAE sign conventions, Fz is negative
fnval(Mz_from_sa_fz_ia,{6.3, -1200, 1.4})

%% Generate higher resolution lookup table data
% evaluate the csaps function at the resolution we want in each dimension
% store the data in a structure we can import into a Simulink n-D lookup block
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
sa = sa_min:sa_res:sa_max;
fz = fz_min:fz_res:fz_max;
ia = ia_min:ia_res:ia_max;

[SA, FZ, IA] = ndgrid(sa, fz, ia);

mz = zeros(size(SA));  % Preallocate the output matrix

% Loop through each combination of inputs and evaluate the function
for i = 1:numel(SA)
    mz(i) = fnval(Mz_from_sa_fz_ia, {SA(i), FZ(i), IA(i)});
end

% validate lookup table data before saving
% plots the same surface we evaluated the csaps spline for above
ia_value = 0;
[~, ia_index] = min(abs(ia - ia_value));  % Find the index of the closest 'ia_value'
mz_fixed_ia = mz(:, :, ia_index)';  % Extract 2D matrix for fixed 'ia_value', and transpose matrix to match dimensions
[SA_grid, FZ_grid] = meshgrid(sa, fz);

% Using surf for a 3D surface plot
figure;
surf(SA_grid, FZ_grid, mz_fixed_ia);
xlabel('sa');
ylabel('fz');
zlabel('mz');
title(['Hi Res Mz Lookup for ia = ', num2str(ia_value)]);

% Validate by comparing values at a point:
fprintf("Value of Csaps: %d\n", fnval(Mz_from_sa_fz_ia, {6, -1600, 0}));
fprintf("Value of Lookup: %d\n", mz(77, 1, 1));

% save both spline and lookup data in .mat file for use in simulink
save("../tire_data/d2704/d2704_7in_csaps_mz.mat", 'Mz_from_sa_fz_ia', 'sa', 'fz', 'ia', 'mz');

% Note: need to debug why plots are not plotting the same (data is fine)
figure('Name','Overturning Moment (Mz) vs. SA, Fz, IA; surf = csaps, dots = hi res lookup')
fnplt(Mz_from_sa_fz_ia)
hold on;
surf(SA_grid, FZ_grid, mz_fixed_ia);
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force (N)')
view(45,45)
%}