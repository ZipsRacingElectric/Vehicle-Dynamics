%{

%% Overview:
This script generates takes in a single slip angle sweep of tire data and
perfoms data cleaning / filtering as an example and to develop the
algorithms.

%% Notes:
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- SAE sign convention 
- sampled at 100 Hz

%% TODO:
- eliminate hysterisis in slip angle sweeps
- fix CSAPS accuracy error (apparent in Mz final plot)

%% Data and units:
AMBTMP: Ambient room temp, deg C
ET: elapsed time, seconds
FX: longituidinal force, N
FY: lateral force, N
FZ: normal load, N
IA: inclination angle, degrees
MX: overturning moment, N-m
MZ: aligning torque, N-m
N: rotational speed, RPM
NFX: normalized longitudinal force (FX/FZ), unitless
NFY: normalized lateral force (FY/FZ), unitless
P: tire pressure, kPa (converted to psi in this script)
RE: effective radius, cm
RL: loaded radius, cm
RST: road surface temp, deg C
SA: slip angle, degrees
SL: slip ratio based on RE (textbook slip ratio), unitless
SR: slip ratio based on RL (used for calspan control), unitless
TSTC: tire surface temp center, deg C
TSTI: tire surface temp inboard, deg C
TSTO: tire surface temp outboard, deg C
V: road speed, kph

%}

clear; close all; clc;
test_data = load("../tire_data/d2704/pure_cornering_test_sweep.mat");

fy = test_data.fy;
fz = test_data.fz;
ia = test_data.ia;
mx = test_data.mx;
mz = test_data.mz;
p = test_data.pressure;
rl = test_data.rl;
sa = test_data.sa;
v = test_data.v;

z = linspace(0, (length(fy)-1), length(fy))'; % create an array of sampled points for plotting
time_step = 1/100; % 100 Hz time sampling
t = z * time_step; % create an array of time starting at the beginning of the sampled data

%% Filtering the Fy Channel
% for each transient sweep, Fz is supposed to be held constant, so we scale Fy
% by multiplying it by the ratio between the Fz mean and the instantaneous Fz.
% since this entire sweep is going to be data for only one Fz point in our
% tire model, this is important

%% Note: filtering this way was found not to be all that helpful

%{
filtering_ratio = mean(fz) ./ fz;
fy_filtered = fy .* filtering_ratio; % element-wise divide


% Plot the filtered data
figure;
hold on;
plot(t, fy);
plot(t, fy_filtered);
plot(t, fz);
hold off;
title("Fy Filtering for Fz Noise");
xlabel("Time (s)");
ylabel("Fy (N)");
xlim([0, t(end)]);
legend("Fy Raw", "Fy Normalized", "Fz");

fy = fy_filtered;
%}


%% Filtering MZ
% Next step is to capture the rational data between the max and minimum
% values, peek at the endpoints, fix up some problems at the MZ endpoints,
% and proceed with data fitting.

mz_filtered = mz;
threshold = 7; % threshold for identifying outliers

[~,imn]=min(sa); % capture slip angle minimum point
[~,imx]=max(sa); % capture slip angle maximum point
p=1:length(sa);
rng=imx-50:imx+50; % This is a range of our data at maximum MZ (used for MZ filtering)

% Being careful not to use a Matlab reserved word for a variable name.
warning off % Keep down the chatter over multiple observations

% fit this data to a polynomial. Crude but fair. We are only using it
% to look for outliers.
pp=polyfit(p(rng),mz(rng)',3);
warning on
mzf = polyval(pp,p(rng)); % filtered mz is the polynomial created

% This step spots data values for MZ that are greater than an arbitray
% level. I believe these spikes are related to the MZ transient response.
% A smarter approach would be to use normalized residuals, but who did
% I just hear volunteer for that task?
ind=find(abs(mzf-mz(rng)') > threshold); % find indexes where the deviation of the filtered data and mz is greater than 7
mz_filtered(rng(ind))=mzf(ind); % set the mz data to the filtered data at these points

% do the same for the minimum Mz point
rng=imn-50:imn+50;% This is a range of our data at minimum MZ
warning off
pp=polyfit(p(rng),mz(rng)',3);
warning on
mzf=polyval(pp,p(rng));
ind=find(abs(mzf-mz(rng)') > threshold);
mz_filtered(rng(ind))=mzf(ind);

% Plot the filtered data
figure;
hold on;
plot(t, mz);
plot(t, mz_filtered);
title("Mz Filtering");
xlabel("Time (s)");
ylabel("Mz (N-m)");
xlim([0, t(end)]);
hold off;
legend("Fy Raw", "Fy Normalized");

%% Remove hysterisis from SA sweep and calculate a 1st order time constant
% we want to compare the hysterisis in Fy values between the increasing
% SA sweep and the decreasing SA sweep
% each SA sweep does this: 0 deg -> +12 deg -> -12 deg -> 0 deg

% organize the data into two arrays representing increasing and
% decreasing slip angle sweeps
[~, n_at_sa_max] = max(sa); % find index of maximum sa
[~, n_at_sa_min] = min(sa); % find index of minimum sa

fy_increasing = cat(1, fy(n_at_sa_min:end), fy(1:n_at_sa_max));
sa_increasing = cat(1, sa(n_at_sa_min:end), sa(1:n_at_sa_max));
fy_decreasing = flip(fy(n_at_sa_max:n_at_sa_min)); % Fy must be "moving" in the same direction as the other array for xcorr to work
sa_decreasing = flip(sa(n_at_sa_max:n_at_sa_min));

% Define a common SA grid
sa_common = linspace(floor(min(sa)), ceil(max(sa)), 100); % Common SA grid

% ensure all sa points are unique before interpolating
[sa_increasing, ~, idx_inc] = unique(sa_increasing);
fy_increasing = accumarray(idx_inc, fy_increasing, [], @mean);
[sa_decreasing, ~, idx_dec] = unique(sa_decreasing);
fy_decreasing = accumarray(idx_dec, fy_decreasing, [], @mean);

% Interpolate Fy onto the common SA grid
Fy_inc_interp = interp1(sa_increasing, fy_increasing, sa_common, 'spline', 'extrap');
Fy_dec_interp = interp1(sa_decreasing, fy_decreasing, sa_common, 'spline', 'extrap');

% Plot results to verify alignment
figure;
plot(sa_increasing, fy_increasing, 'r-', sa_decreasing, fy_decreasing, 'b--');
hold on;
plot(sa_common, Fy_inc_interp, 'ro', sa_common, Fy_dec_interp, 'bo');
legend('Fy Increasing (Original)', 'Fy Decreasing (Original)', ...
       'Fy Increasing (Interpolated)', 'Fy Decreasing (Interpolated)');
xlabel('Slip Angle (SA)');
ylabel('Lateral Force (Fy)');
title('Interpolated Fy Data');
hold off;

% plot fy vs sa
figure;
plot(sa, fy)

% determine the time constant for relaxation length
[cross_corr, lags] = xcorr(fy_increasing, fy_decreasing, 'none');
[~, max_index] = max(cross_corr);
time_lag = lags(max_index) * time_step; % this is approximately the time constant for 1st order relaxation length

% plot cross correlation stuff
figure;
plot(lags, cross_corr);
xlabel("Lag");
ylabel("Correlation");
title("Cross correlation for Fy, increasing and decreasing SA sweeps")

%% Frequency Analysis of data
% first plot the FFT of the important channels to determine high frequency
% noise

mz_spectrum = fft(mz);
fy_spectrum = fft(fy);

L = length(mz);
freq = 1/time_step; % sample rate

figure;
plot(freq / L * (0:L-1), abs(mz_spectrum))
title("Spectrum of Mz channel")
xlabel("f (Hz)");
ylabel("|fft(Mz)|")
xlim([0 10]);

figure;
plot(freq / L * (0:L-1), abs(fy_spectrum))
title("Spectrum of Fy channel")
xlabel("f (Hz)");
ylabel("|fft(Fy)|")
xlim([0 10]);

%% noise clearly identified in Mz at 4 & 7 Hz

%% Low Pass Filter on Mz
% 3rd order finite impulse response filter
filtertype = 'FIR';
N = 3;
freq_pass = 2; % attentuate frequencies greater than 3 Hz
freq_stop = 10;
Rp = 0.5;
Astop = 50;

LPF_mz = dsp.LowpassFilter('SampleRate', freq, 'PassbandFrequency', freq_pass, 'StopbandFrequency', freq_stop, 'PassbandRipple', Rp, 'StopbandAttenuation', Astop);
mz_filtered = step(LPF_mz, mz);

% with the phase delay introduced by the filter, get the cross-correlation
% and shift back the data the correct amount
[cross_corr, lags] = xcorr(mz, mz_filtered, 'none');
[~, max_index] = max(cross_corr);
sample_lead = -lags(max_index);

% shift the data back to where it belongs
mz_filtered_shifted = mz_filtered(sample_lead:end);

% plot cross correlation
figure;
plot(lags, cross_corr);
xlabel("Lag");
ylabel("Correlation");
title("Cross correlation for Mz and Mz filtered")

% Plot the filtered data
figure;
hold on;
plot(t, mz);
plot(t, mz_filtered);
plot(t(1:length(mz_filtered_shifted)), mz_filtered_shifted);
title("Mz Filtering via Low Pass Filter");
xlabel("Time (s)");
ylabel("Mz (N-m)");
xlim([0, t(end)]);
hold off;
legend("Mz Raw", "Mz Filtered", "Mz Filtered and Shifted");

%% Spline fitting the continuous data to subset it with 1 degree slip angle increments
% with some tighter tension:
sp_fy=csaps(sa,fy,.1);
sp_mz=csaps(sa,mz,.1);
sp_mx=csaps(sa,mx,.1);
sp_rl=csaps(sa,rl,.1);