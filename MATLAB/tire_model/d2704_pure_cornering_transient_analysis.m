%{
%% Overview:
This script takes a look at the transient cornering test plans to help 
determing relaxation length effects. 
It has been developed from Bill Cobb's information (availible on TTC forum)

For use with round 9 data.

%% Notes:
- You will need the tire_data folder from the OneDrive. Put it in the
MATLAB directory of the local Vehicle-Dynamics-ZR25 repository
- For use with round 9 data.
- SAE sign convention 
- sampled at 100 Hz
- use .MAT files with SI units for this script
- each TTC round uses a different test procedure, so this script needs
modified for anything that isn't round 9

%% TODO:
- add relaxation length estimates from slip angle sweep hysterisis to data
set

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
%{
The Round 9 transient test runs are identical to those in Rounds 6-8. They 
begin with initial spring rate tests. Then, several low-speed transient 
tests are performed, all at zero inclination angle. See Figure 1. Note that
 all subplots in Figure 1 have the same x-axis time scale.

Transient tests at each operating condition (load and inflation pressure) 
involve stopping the tire, steering to a new slip angle at zero speed and 
then rolling the tire forward slowly. Steer angles of −1, +1 and +6 deg are
 used, with a return to 0 deg for each step steer. As the tire rolls slowly
forward, the force and moment channels respond as the tire seeks a new 
equilibrium condition, see Figure 2. This response distance, sometimes 
characterized by a “relaxation length” metric, is an indication of tire 
transient response.
%}

load("../tire_data/d2704/data/B2356run22.mat");

%% Clean tire data
% delete initial spring rate tests, only care about relaxation length rn
p = 1375; % data point which we want to cut off to
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

% convert tire pressure data from kPa to PSI
P = P * 0.14503773773020923;

% convert velocity from kph to m/s
V = V * (1000/1) * (1/60 * 1/60);


%% Plot raw data
figure;

subplot(6,1,1)
plot(ET, P)
title('Pressure (psi)')
xlabel('Elapsed Time (s)')
ylabel('Pressure (psi)')

subplot(6,1,2)
plot(ET, SA)
title("Slip Angle (deg)")
xlabel('Elapsed Time (s)')
ylabel('Slip Angle (deg)')

subplot(6,1,3)
plot(ET, FZ)
title("Normal Load (N)")
xlabel('Elapsed Time (s)')
ylabel('Normal Load (N)')

subplot(6,1,4)
plot(ET, FY)
title("Lateral Force (N)")
xlabel('Elapsed Time (s)')
ylabel('Lateral Force (N)')

subplot(6,1,5)
plot(ET, IA)
title("Inclination Angle (deg)")
xlabel('Elapsed Time (s)')
ylabel('Inclination Angle (deg)')

subplot(6,1,6)
plot(ET, V)
title("Velocity (m/s)")
xlabel('Elapsed Time (s)')
ylabel('Velocity (m/s)')

%% Locate start / end of each individual transient sweep
% create a spline equation to locate zero crossing points in data velocity
% data, which represent the start of each transient test
m = 1:length(V);

% V is subtracted by 0.3 m/s so zero crossings aren't identified in the 0 m/s noise
sp = spline(m,V - 0.3);

% location of zero crossings in data
z = fnzeros(sp);
z = round(z(1,:));

% plot raw slip angle data and zero crossing points
figure;
hold on
plot(m, V)

% plot the spline and zeros
plot(z,zeros(length(z)),'bo')
xlabel('Point Count')
ylabel('Velocity')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff

% Don't want to remember any data from apprevious processing session:
clear fmdata;


q = 0;
% Subset the data into individual transient sweeps:

for n=1:2:(length(z)-2) % every 3 zeros represents a transient sweep

    %% capture the individual channel data between the zero crossings.
    % the data is only valid while velocity is moving, so grab data up
    % until velocity drops back to zero
    et=ET(z(n):z(n+1));
    v=V(z(n):z(n+1));
    sa=SA(z(n):z(n+1)); 
    fz=FZ(z(n):z(n+1));
    fy=FY(z(n):z(n+1));
    mz=MZ(z(n):z(n+1));
    mx=MX(z(n):z(n+1));
    ia=IA(z(n):z(n+1));
    pressure=P(z(n):z(n+1));

    %% normalize the Fy channel to account for Fz noise
    % for each transient sweep, Fz is held constant, so we take the
    % average:
    fy = (fy * mean(fz)) ./ fz; % element-wise divide

    %% fit a 1st order curve to data
    % get initial value
    fy_i = fy(1);
    
    % get steady state value by averaging the latter half of the data
    % (assumption: transient is dissipated already in this portion of fy)
    fy_ss = mean(fy(ceil(length(fy)/2)+1:end));

    % shift the time array so it starts at t = 0 for first order fitting
    et_shift = et - et(1);

    % Define the first-order system model as an anonymous function of tau
    % and t
    first_order_model_rise = @(tau, t) (1 - exp(-t / tau));
    first_order_model_decay = @(tau, t) (exp(-t / tau));

    % Initial guesses for tau
    tau_init = et_shift(end) / 10; % Initial guess for time constant

    % shift data then scale between 0 and 1
    if fy_ss < fy_i % if we are steering out for this transient

        fy_norm = (fy - fy_ss) / (fy_i - fy_ss);

        % Nonlinear curve fitting  on normalized fy to get optimal tau
        tau = lsqcurvefit(first_order_model_decay, tau_init, et_shift, fy_norm);
        first_order_fy = first_order_model_decay(tau, et_shift);

    else % if we are steering in for this transient
        fy_norm = (fy - fy_i) / (fy_ss - fy_i);

        % Nonlinear curve fitting  on normalized fy to get optimal tau
        tau = lsqcurvefit(first_order_model_rise, tau_init, et_shift, fy_norm);
        first_order_fy = first_order_model_rise(tau, et_shift);

    end

    %% Spline fit the continuous data
    sp_v=csaps(et,v,.95);
    sp_sa=csaps(et,sa,.95);
    sp_fy=csaps(et,fy_norm,.8);
    sp_mz=csaps(et,mz,.95);
    sp_mx=csaps(et,mx,.95);

    %% Check out Segment 9
    % Just out of curiosity, what kind of data are we dealing with?
    if isequal(n,19)


        fprintf(num2str(tau))

        figure;
        subplot(4,1,1)
        hold on
        plot(et,fy,'.','color',[.5 .5 .5])
        title({['Fz= ' num2str(round(mean(fz))) ' N'];['IA= ' num2str(round(mean(ia))) '°, ' num2str(round(mean(pressure))) ' psi']})
        xlabel('Elapsed Time, sec')
        ylabel('Norm Lateral Force Fy/Fz')
        legend({['Test Data','tau = ' num2str(tau)]})

        subplot(4,1,2)
        hold on
        plot(et,fz,'.','color',[.5 .5 .5])
        xlabel('Elapsed Time, sec')
        ylabel('Normal Force Fz, N')
        
        subplot(4,1,3)
        hold on
        plot(et,sa,'.','color',[.5 .5 .5])
        fnplt(sp_sa,'b')
        xlabel('Elapsed Time, sec')
        ylabel('Slip Angle, deg')

        subplot(4,1,4)
        hold on
        plot(et,v,'.','color',[.5 .5 .5])
        fnplt(sp_v,'b')
        xlabel('Elapsed Time, sec')
        ylabel('Velocity, m/s') 
        hold off

        % Plot the original data and fitted curve
        figure;
        plot(et_shift, fy_norm, 'bo', 'DisplayName', 'Original Data'); hold on;
        plot(et_shift, first_order_fy, 'r-', 'DisplayName', ['Fitted Curve (Tau = ', num2str(tau), ')']);
        xlabel('Elapsed Time (s)');
        ylabel('Normal Force (Fy)');
        title('First-Order System Fitting');
        legend show;
        grid on;

    end
    
    
    %% create a clean, fitted dataset called fmdata
    % columns: SA, IA, Fz, Fy, Mz, Mx, P
    for sl=floor(min(sa)):1:ceil(max(sa)); % full sweep of slip angle spline in whole integers of degrees
        q=q+1;
        fmdata(q,1)=sl; 
        fmdata(q,2)=round(mean(ia));
        fmdata(q,3)=mean(fz);
        fmdata(q,4)=fnval(sp_fy,sl);
        fmdata(q,5)=fnval(sp_mz,sl);
        fmdata(q,6)=fnval(sp_mx,sl);
        fmdata(q,7)=round(mean(pressure));
    
    end 
    %}
end

%{
%% Sort Data Array 1st by Pressure, Camber, Slip, and finally Fz
% Use |sortrows| to preserve the array correspondence.
fmdata = sortrows(fmdata,[7,2,1,3]);
incls = unique(round(fmdata(:,2)))';
nincls = length(incls); % number of distinct inclination angles
slips = unique(round(fmdata(:,1)))';
nslips = length(slips); % number of distinct slip angles
press = unique(round(fmdata(:,7)))';
npress = length(press); % number of distinct pressures

%% Save organized tire data
save("d2704_7in_fmdata.mat", "fmdata");

%% Data Sets (Slip, Load, Camber, Pressure):
inx0 = find(fmdata(:,2) == 0); % array indexes with zero camber points
fmdata0 = fmdata(inx0,:); % selects rest of data corresponding to 0 IA
px0 = find(fmdata0(:,7) == 8); % array indexes for zero camber and 8 psi
fmdata0ip = fmdata0(px0,:); % selects rest of 0 IA data corresponding to 8 psi

% Next, we transpose the arrays to make our spline functions happy:
reshaped_fmdata = reshape(fmdata0ip(:,3),[],nslips);
loads = mean(reshape(fmdata0ip(:,3),[],nslips),2)'; % appears to just give us a list of unique loads
nloads = length(loads);

% Take a look at FZ:
fz0ip = reshape(fmdata0ip(:,3),nloads,nslips)'; % create an array with nload rows x nslip columns
fy0ip = reshape(fmdata0ip(:,4),nloads,nslips)'
mz0ip = reshape(fmdata0ip(:,5),nloads,nslips)';
mx0ip = reshape(fmdata0ip(:,6),nloads,nslips)';

% normalized lateral force (mu value)
nfy0ip = fy0ip./fz0ip;

%{
%% FY Surface Fit (0 IA, 8 psi)
LATE_SLIP_VERT = csaps({slips,loads},fy0ip)
figure('Name','Lateral Force vs. Slip Angle & Vertical Load')
fnplt(LATE_SLIP_VERT)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force (N)')
view(45,45)

CS=fnder(LATE_SLIP_VERT,[1,0])
figure('Name','Cornering Stiffness vs. Slip Angle & Vertical Load')
fnplt(CS)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Cornering Stiffness (N/deg)')

%% Normalized Lateral Force Surface Fit (0 IA, 8 psi)
NLATE_SLIP_VERT = csaps({slips,loads},nfy0ip)
figure('Name','Load Normalized Lateral Force vs. Slip Angle & Vertical Load ')
fnplt(NLATE_SLIP_VERT)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Lateral Force mu value (Fy/Fz)')
view(45,45)
% Wow, look at the mu on that baby, good as a Sprint Cup Left Side Tire !!
% Here's the traditional normalized cornering stiffness used in industry: 
NCS=fnder(NLATE_SLIP_VERT,[1,0])
figure('Name', 'Normalized Cornering Stiffness vs. Slip Angle & Vertical Load ')
fnplt(NCS)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Normalized Cornering Stiffness (N/deg/N)')

%% Mz Surface Fit (0 IA, 8 psi)
ALNT_SLIP_VERT = csaps({slips,loads},mz0ip)
figure('Name','Aligning Moment vs. Slip Angle & Vertical Load')
fnplt(ALNT_SLIP_VERT)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Aligning Moment Mz (Nm/deg)')

%% Mx Surface Fit (0 IA, 8 psi)
OVTM_SLIP_VERT = csaps({slips,loads},mx0ip)
figure('Name','Overturning Moment vs. Slip Angle & Vertical Load')
fnplt(OVTM_SLIP_VERT)
view(30,45)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Overturning Moment Mx (Nm)')

%% Pneumatic Scrub
PSCRUB_SLIP_VERT = csaps({slips,loads},1000*mx0ip./fz0ip) 
figure('Name','Pneumatic Scrub vs. Slip Angle & Vertical Load')
fnplt(PSCRUB_SLIP_VERT)
view(30,45)
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Pneumatic Scrub (mm)')

%% Pneumatic Trail Surface
PTRAIL_SLIP_VERT = csaps({slips,loads},1000*mz0ip./fy0ip,.707)
figure('Name','Pneumatic Trail vs. Slip Angle & Vertical Load')
fnplt(PTRAIL_SLIP_VERT)
view(30,45)
title('Although I would not produce Pneumatic Trail this way, good guess, though ...' )
xlabel('Slip Angle (deg)')
ylabel('Vertical Load (N)')
zlabel('Pneumatic Trail (mm)')

%% for graphs on Fy = 0 surfaces, see Matlab Processing of FSAE TTC Tire Test Data.pdf
%}
%}