%{

%% Overview:
This script generates formatted data sets from the raw TTC files. This data
can then be used with other scripts to generate tire models from. 
It has been modified from Bill Cobb's matlab script 
(availible on TTC forum) for use with round 9 data.

%% Notes:
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- When you run this script, the formatted data is saved in the tire directory as d2704_7in_formatted_data.mat
- SAE sign convention 
- sampled at 100 Hz
- use .MAT files with SI units for this script
- each TTC round uses a different test procedure, so this script needs
modified for anything that isn't round 9

%% TODO:
- Filter, interpolate, & shift Mz using FIR developed in data_cleaning_example.m
- eliminate hysterisis in slip angle sweeps

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

%% Load tire data
% pure slip 12,10, and 14 psi. This is run data so the
% break in procedures, spring rates, have already been removed.
% Note: the 12 psi sweep in this file is a new tire, the second file
% contains a 12 psi sweep at a more broken in tire
file1 = load("../tire_data/d2704/data/B2356run23.mat");
file2 = load("../tire_data/d2704/data/B2356run24.mat");

% Load low speed relaxation length data
%% load pure cornering fitted data points
load("../tire_data/d2704/d2704_7in_relaxation_low_speed.mat");

time_step = 1 / 100; % time step between data points

% concatenate the data into single channels
AMBTMP = [file1.AMBTMP; file2.AMBTMP];
ET = [file1.ET; (file2.ET + file1.ET(length(file1.ET)))]; % shift ET time by the length of the first file
FX = [file1.FX; file2.FX];
FY = [file1.FY; file2.FY];
FZ = [file1.FZ; file2.FZ];
IA = [file1.IA; file2.IA];
MX = [file1.MX; file2.MX];
MZ = [file1.MZ; file2.MZ];
N = [file1.N; file2.N];
NFX = [file1.NFX; file2.NFX];
NFY = [file1.NFY; file2.NFY];
P = [file1.P; file2.P];
RE = [file1.RE; file2.RE];
RL = [file1.RL; file2.RL];
RST = [file1.RST; file2.RST];
RUN = [file1.RUN; file2.RUN];
SA = [file1.SA; file2.SA];
SL = [file1.SL; file2.SL];
SR = [file1.SR; file2.SR];
TSTC = [file1.TSTC; file2.TSTC];
TSTI = [file1.TSTI; file2.TSTI];
TSTO = [file1.TSTO; file2.TSTO];
V = [file1.V; file2.V];

%plot SA vs ET to verify the data is good before procceding
%figure;
%plot(ET, P)

%% Clean tire data
% delete the first 12 psi run, because it is a new tire and the 2nd 12 psi
% run should be more accurate (broken in and at temp)
p = 21170; % data point which we want to cut off to
AMBTMP(1:p,:)=[];
ET(1:p,:)=[];
FX(1:p,:)=[];
FY(1:p,:)=[];
FZ(1:p,:)=[];
IA(1:p,:)=[];
MX(1:p,:)=[];
MZ(1:p,:)=[];
N(1:p,:)=[];
NFX(1:p,:)=[];
NFY(1:p,:)=[];
P(1:p,:)=[];
RE(1:p,:)=[];
RL(1:p,:)=[];
RST(1:p,:)=[];
RUN(1:p,:)=[];
SA(1:p,:)=[];
SL(1:p,:)=[];
SR(1:p,:)=[];
TSTC(1:p,:)=[];
TSTI(1:p,:)=[];
TSTO(1:p,:)=[];
V(1:p,:)=[];

% 898.25 s
% get rid of the first SA sweep for 10 psi. this appears to be a duplicate
% 0 IA, -1100 N condition
p = 1253; % data point which we want to cut off to
AMBTMP(1:p,:)=[];
ET(1:p,:)=[];
FX(1:p,:)=[];
FY(1:p,:)=[];
FZ(1:p,:)=[];
IA(1:p,:)=[];
MX(1:p,:)=[];
MZ(1:p,:)=[];
N(1:p,:)=[];
NFX(1:p,:)=[];
NFY(1:p,:)=[];
P(1:p,:)=[];
RE(1:p,:)=[];
RL(1:p,:)=[];
RST(1:p,:)=[];
RUN(1:p,:)=[];
SA(1:p,:)=[];
SL(1:p,:)=[];
SR(1:p,:)=[];
TSTC(1:p,:)=[];
TSTI(1:p,:)=[];
TSTO(1:p,:)=[];
V(1:p,:)=[];

% 1492 to 1505
% get rid of extra 1100 N condition
p1 = 18697; % data point which we want to start cut off to
p2 = 19942; % data point which we want to cut off to
AMBTMP(p1:p2,:)=[];
ET(p1:p2,:)=[];
FX(p1:p2,:)=[];
FY(p1:p2,:)=[];
FZ(p1:p2,:)=[];
IA(p1:p2,:)=[];
MX(p1:p2,:)=[];
MZ(p1:p2,:)=[];
N(p1:p2,:)=[];
NFX(p1:p2,:)=[];
NFY(p1:p2,:)=[];
P(p1:p2,:)=[];
RE(p1:p2,:)=[];
RL(p1:p2,:)=[];
RST(p1:p2,:)=[];
RUN(p1:p2,:)=[];
SA(p1:p2,:)=[];
SL(p1:p2,:)=[];
SR(p1:p2,:)=[];
TSTC(p1:p2,:)=[];
TSTI(p1:p2,:)=[];
TSTO(p1:p2,:)=[];
V(p1:p2,:)=[];

% 2121 to 2133
% get rid of extra 1100 N condition
p1 = 37379; % data point which we want to start cut off to
p2 = 38625; % data point which we want to cut off to
AMBTMP(p1:p2,:)=[];
ET(p1:p2,:)=[];
FX(p1:p2,:)=[];
FY(p1:p2,:)=[];
FZ(p1:p2,:)=[];
IA(p1:p2,:)=[];
MX(p1:p2,:)=[];
MZ(p1:p2,:)=[];
N(p1:p2,:)=[];
NFX(p1:p2,:)=[];
NFY(p1:p2,:)=[];
P(p1:p2,:)=[];
RE(p1:p2,:)=[];
RL(p1:p2,:)=[];
RST(p1:p2,:)=[];
RUN(p1:p2,:)=[];
SA(p1:p2,:)=[];
SL(p1:p2,:)=[];
SR(p1:p2,:)=[];
TSTC(p1:p2,:)=[];
TSTI(p1:p2,:)=[];
TSTO(p1:p2,:)=[];
V(p1:p2,:)=[];

% 2733 to 2745
% get rid of extra 1100 N condition
p1 = 56066; % data point which we want to start cut off to
p2 = 57310; % data point which we want to cut off to
AMBTMP(p1:p2,:)=[];
ET(p1:p2,:)=[];
FX(p1:p2,:)=[];
FY(p1:p2,:)=[];
FZ(p1:p2,:)=[];
IA(p1:p2,:)=[];
MX(p1:p2,:)=[];
MZ(p1:p2,:)=[];
N(p1:p2,:)=[];
NFX(p1:p2,:)=[];
NFY(p1:p2,:)=[];
P(p1:p2,:)=[];
RE(p1:p2,:)=[];
RL(p1:p2,:)=[];
RST(p1:p2,:)=[];
RUN(p1:p2,:)=[];
SA(p1:p2,:)=[];
SL(p1:p2,:)=[];
SR(p1:p2,:)=[];
TSTC(p1:p2,:)=[];
TSTI(p1:p2,:)=[];
TSTO(p1:p2,:)=[];
V(p1:p2,:)=[];

% plot SA vs ET to verify the data is good before procceding
%{
figure;
plot(ET, P)
figure;
plot(ET, SA)
%}

% convert tire pressure data from kPa to PSI
P = P * 0.14503773773020923;

% plot SA vs ET to verify the data is good before procceding
% important: made sure the data crosses zero before and after the SA sweeps
%{
figure;
plot(ET, P)
figure;
plot(ET, SA)
figure;
hold on;
yyaxis left;
plot(ET, SA)
yyaxis right;
plot(ET, FZ)
hold off;
figure;
plot(ET, IA)
%}

%% Subset the data into individual SA sweeps:
% create a spline equation to locate zero crossing points in data
m = 1:length(SA);
sp = spline(m,SA);

% location of zero crossings in data
z = fnzeros(sp);
z = round(z(1,:));

% plot raw slip angle data and zero crossing points
figure;
hold on
plot(m, SA)

% plot the spline and zeros
plot(z,zeros(length(z)),'bo')
xlabel('Point Count')
ylabel('Slip Angle')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff

q = 0;
for n=1:2:(length(z)-2) % every 3 zero crossings represents a full SA sweep

    %% capture the individual channel data:
    sa=SA(z(n):z(n+2)); % n to n+2 represents going through a full (+-) 12 deg SA sweep
    fz=FZ(z(n):z(n+2));
    fy=FY(z(n):z(n+2));
    mz=MZ(z(n):z(n+2));
    mx=MX(z(n):z(n+2));
    rl=RL(z(n):z(n+2));
    ia=IA(z(n):z(n+2));
    pressure=P(z(n):z(n+2));
    v=V(z(n):z(n+2));
    % Now we have collected the tire channels for each full slip sweep.

    %% filter the Fy channel to account for Fz noise
    % for each transient sweep, Fz is held constant, so we filter Fy for Fz
    % noise contamination by multiplying it by the ratio between the Fz mean and the instantaneous Fz.
    % fy = (fy * mean(fz)) ./ fz; % element-wise divide
    %% Note: for certain sweeps this was not helpful, disabled for now
    
    %% Filter MZ data
    % Next step is to capture the rational data between the max and minimum
    % values, peek at the endpoints, fix up some problems at the MZ endpoints,
    % and proceed with data fitting.

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
    mzf=polyval(pp,p(rng));

    % This step spots data values for MZ that are greater than an arbitray
    % level. I believe these spikes are related to the MZ transient response.
    % A smarter approach would be to use normalized residuals, but who did
    % I just hear volunteer for that task?
    ind=find(abs(mzf-mz(rng)') > 7);
    mz(rng(ind))=mzf(ind);
    rng=imn-50:imn+50;% This is a range of our data at minimum MZ
    warning off
    pp=polyfit(p(rng),mz(rng)',3);
    warning on
    mzf=polyval(pp,p(rng));
    ind=find(abs(mzf-mz(rng)') > 7);
    mz(rng(ind))=mzf(ind);

    %% Calculate 1st order relaxation length time constant
    % we want to compare the hysterisis in Fy values between the increasing
    % SA sweep and the decreasing SA sweep
    % each SA sweep does this: 0 deg -> +12 deg -> -12 deg -> 0 deg

    % organize the data into two arrays representing increasing and
    % decreasing slip angle sweeps
    [~, n_at_sa_max] = max(sa); % find index of maximum sa
    [~, n_at_sa_min] = min(sa); % find index of minimum sa

    fy_increasing_sa = cat(1, fy(n_at_sa_min:end), fy(1:n_at_sa_max));
    fy_decreasing_sa = flip(fy(n_at_sa_max:n_at_sa_min)); % Fy must be "moving" in the same direction as the other array for xcorr to work

    [cross_corr, lags] = xcorr(fy_increasing_sa, fy_decreasing_sa, 'none');
    [~, max_index] = max(cross_corr); % find index for best correlation
    tau = lags(max_index) * time_step; % time lag associated with best correlation
    % this time lag is a good approximation for 1st order relaxation length
    % time constant
    % Possibly better: fit a 1st order equation to the data using
    % lsqcurvefit

    %% Spline fitting the continuous data to subset it with 1 Degree slip angle increments
    % with some tighter tension:
    sp_fy=csaps(sa,fy,.1);
    sp_mz=csaps(sa,mz,.1);
    sp_mx=csaps(sa,mx,.1);
    sp_rl=csaps(sa,rl,.1);

    %% Check out Segment 17
    % Just out of curiosity, what kind of data are we dealing with?
    if isequal(n,17)

        % save raw data
        save("../tire_data/d2704/pure_cornering_test_sweep.mat", "sa", "fz", "fy", "mz", "mx", "rl", "ia", "pressure", "v");

        fprintf("hit segment 9")
        figure;
        subplot(3,1,1)
        hold on
        plot(sa,fy,'.','color',[.5 .5 .5])
        fnplt(sp_fy,'b')
        title({['Fz= ' num2str(round(mean(fz))) ' N'];['IA= ' num2str(round(mean(ia))) '°, ' num2str(round(mean(pressure))) ' psi']})
        xlabel('Slip Angle')
        ylabel('Lateral Force')
        line([min(sa) max(sa)],[0 0],'color','k')
        line([0 0],[min(fy) max(fy)],'color','k')
        legend('Test Data','Fitted Data')
        
        subplot(3,1,2)
        hold on
        plot(sa,mz,'.','color',[.5 .5 .5])
        fnplt(sp_mz,'b')
        xlabel('Slip Angle')
        ylabel('Aligning Moment')
        line([min(sa) max(sa)],[0 0],'color','k')
        line([0 0],[min(mz) max(mz)],'color','k')
        subplot(3,1,3)
        hold on
        plot(sa,mx,'.','color',[.5 .5 .5])
        fnplt(sp_mx,'b')
        xlabel('Slip Angle')
        ylabel('Overturning Moment') 
        line([min(sa) max(sa)],[0 0],'color','k') 
        line([0 0],[min(mz) max(mz)],'color','k')

        % plot cross correlation stuff
        figure;
        plot(lags, cross_corr);
        xlabel("Lag");
        ylabel("Correlation");
        title("Cross correlation for Fy, increasing and decreasing SA sweeps")
    end

    %% create formatted data
    % Each column: time constant, velocity, normal force, pressure
    q = q + 1;
    relaxation_high_speed(q,1) = tau; % for fy only
    relaxation_high_speed(q,2) = round(mean(v));
    relaxation_high_speed(q,3) = mean(fz);

end

%% Combine Low Speed data and format

% eliminate outliers in low speed data
outlier_indexes = cat(1, find(relaxation_low_speed(:, 1) > 2), find(relaxation_low_speed(:, 1) < 0));
relaxation_low_speed(outlier_indexes, :) = []; % remove row

% Clean up low speed data
% for simplicity assume equal relaxation coefficient across all normal
% loads
unique_fz = unique(round(relaxation_low_speed(:, 3)));
low_speed_data = [mean(relaxation_low_speed(:,1)) + zeros(length(unique_fz), 1), 1 + zeros(length(unique_fz), 1), unique_fz];

% combine data
relaxation_data = [low_speed_data; relaxation_high_speed];

% Add a column for calculated relaxation length
relaxation_length = relaxation_data(:, 1) .* relaxation_data(:, 2); % [m]
relaxation_data = [relaxation_data(:, 1), relaxation_length, relaxation_data(:, 2), relaxation_data(:, 3)];

% Add a zero speed condition
% relaxation coefficient at Fz = 0 should be zero to avoid tire force delays
% when the tire is in the air, ensure csaps interpolates correctly
unique_v = unique(round(relaxation_data(:, 3)));
zero_cond = [zeros(length(unique_v), 2), unique_v, zeros(length(unique_v), 1)];

relaxation_data = [relaxation_data; zero_cond];

% Sort Data by velocity, fz
relaxation_data = sortrows(relaxation_data,[3, 4, 1]);

% Save organized tire data
% Columns are: time constant, relaxation length, velocity, normal load
save("../tire_data/d2704/d2704_7in_relaxation_data.mat", "relaxation_data");

%% display relaxation length coefficients
figure;
scatter3(relaxation_data(:, 3), relaxation_data(:, 4), relaxation_data(:, 1));

% Spline fit to coefficients
% Extract columns
relaxation_coeff = relaxation_data(:, 1);
velocity = relaxation_data(:, 3);
normal_force = relaxation_data(:, 4);

% Create a grid for bivariate data
[velocity_grid, force_grid] = ndgrid(unique(velocity), unique(normal_force));

% Map relaxation_coeff onto the grid
relaxation_coeff_grid = griddata(velocity, normal_force, relaxation_coeff, velocity_grid, force_grid, 'linear');

% Fit a csaps bivariate spline for relaxation coefficient
smoothing_param = 0.01; % Adjust smoothing parameter as needed
spline_bivariate = csaps({unique(velocity), unique(normal_force)}, relaxation_coeff_grid, smoothing_param);

% Generate fine grids for plotting
velocity_fine = linspace(min(velocity), max(velocity), 100);
force_fine = linspace(min(normal_force), max(normal_force), 100);
[velocity_fine_grid, force_fine_grid] = ndgrid(velocity_fine, force_fine);

% Evaluate the spline on the fine grid
relaxation_fit = fnval(spline_bivariate, {velocity_fine, force_fine});

% Plot the results
figure;

% 3D Surface Plot
surf(velocity_fine_grid, force_fine_grid, relaxation_fit, 'EdgeColor', 'none');
hold on;
scatter3(velocity, normal_force, relaxation_coeff, 'filled', 'r');
xlabel('Velocity');
ylabel('Normal Force');
zlabel('Relaxation Coefficient');
title('Bivariate Cubic Smoothing Spline');
grid on;
colorbar;
