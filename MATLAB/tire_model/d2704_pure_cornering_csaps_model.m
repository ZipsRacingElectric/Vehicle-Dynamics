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
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- SAE sign convention 
- load the formatted_data.mat file, generated using the data_formatting.m
script in this same directory
- formatted_data expected columns(in order): SA, IA, Fz, Mz, Mx, P

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

%% select data that holds P = 8 psi constant
incls = unique(round(fmdata(:,2)))';
nincls = length(incls); % figure out number of distinct inclination angles
slips = unique(round(fmdata(:,1)))';
nslips = length(slips); % figure out number of distinct slip angles
press = unique(round(fmdata(:,7)))';
npress = length(press); % figure out number of distinct pressures

p8 = find(fmdata(:,7) == 8); % find array indexes for 8 psi
fmdata_p8 = fmdata(p8,:); % selects all data corresponding to 8 psi

%% Next, we transpose the arrays to make our spline functions happy:
% take the sorted Fz column data and reshape it into columns. Each column
% represents a distint slip angle and inclination angle combo
loads = reshape(fmdata_p8(:,3),[],(nslips*nincls));

% take the average of each row to come up with distinct loads. Also
% transpose into a row vector.
loads = mean(loads,2)';
nloads = length(loads);

% reshape the data into 3D load, slip, and IA dimensions representive of the
% surface csaps will try and fit
fy_p8 = reshape(fmdata_p8(:,4),nloads,nslips,nincls);

% permutate (switch) the rows and columns such that we get:
% rows: SA sweep
% colums: Fz sweep
% vertical dimension: IA
fy_p8 = permute(fy_p8, [2 1 3]);

%% FY Surface Fit (P = 8 psi)
Fy_from_sa_fz_ia = csaps({slips,loads,incls},fy_p8);

% note: when plotting the 4D function, it will give you a 3D surface at IA = 0
figure('Name','Lateral Force vs. Slip Angle, Vertical Load, Inclination Angle')
fnplt(Fy_from_sa_fz_ia)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force (N)')
title('Name','Lateral Force vs. Slip Angle, Vertical Load, IA = 0 deg')
view(45,45)

%% evaluating the CSAPS tire model
%keep in mind the SAE sign conventions, Fz is negative
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

fy = zeros(size(SA));  % Preallocate the output matrix

% Loop through each combination of inputs and evaluate the function
for i = 1:numel(SA)
    fy(i) = fnval(Fy_from_sa_fz_ia, {SA(i), FZ(i), IA(i)});
end

%% validate lookup table data before saving
% plots the same surface we evaluated the csaps spline for above
ia_value = 0;
[~, ia_index] = min(abs(ia - ia_value));  % Find the index of the closest 'ia_value'
fy_fixed_ia = fy(:, :, ia_index);  % Extract 2D matrix for fixed 'ia_value'
[SA_grid, FZ_grid] = meshgrid(sa, fz);
fy_fixed_ia = fy_fixed_ia';  % Transpose to match dimensions

% Using surf for a 3D surface plot
figure;
surf(SA_grid, FZ_grid, fy_fixed_ia);
xlabel('sa');
ylabel('fz');
zlabel('fy');
title(['fy for ia = ', num2str(ia_value)]);

%% save both spline and lookup data in .mat file for use in simulink
save("../tire_data/d2704/d2704_7in_csaps.mat", 'Fy_from_sa_fz_ia', 'sa', 'fz', 'ia', 'fy');

figure('Name','Lateral Force vs. Slip Angle, Vertical Load, Inclination Angle')
fnplt(Fy_from_sa_fz_ia)
hold on;
plot3(SA_grid(:), FZ_grid(:), fy_fixed_ia(:), 'o');  % 'o' specifies dots
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force (N)')
view(45,45)