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

%% Load tire data
% combined slip 8,10,14, and 12 psi at hot temp. This is run data so the
% break in procedures, spring rates, have already been removed.
% Note: the first 12 psi sweep in this file is a new tire, the second file
% contains a 12 psi sweep at a more broken in tire. We will remove the
% first 12 psi run
file1 = load("../tire_data/d2704/data/B2356run57.mat");
file2 = load("../tire_data/d2704/data/B2356run58.mat");

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

% convert tire pressure data from kPa to PSI
P = P * 0.14503773773020923;

%plot SA vs ET to verify the data is good before procceding
%{
figure;
plot(ET, P)
%}

%% Clean tire data
% delete the first 12 psi run, because it is a new tire and the 2nd 12 psi
% run should be more accurate (broken in and at temp)
p = 22739; % data point which we want to cut off to
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

% get rid of the final speed test which occurs during the last 12 psi
% tests. A little bit of the speed change is left to ensure the last slip
% ratio peak is captured
p = 92841; % data point which we want to cut off to
AMBTMP(p:end,:)=[];
ET(p:end,:)=[];
FX(p:end,:)=[];
FY(p:end,:)=[];
FZ(p:end,:)=[];
IA(p:end,:)=[];
MX(p:end,:)=[];
MZ(p:end,:)=[];
N(p:end,:)=[];
NFX(p:end,:)=[];
NFY(p:end,:)=[];
P(p:end,:)=[];
RE(p:end,:)=[];
RL(p:end,:)=[];
RST(p:end,:)=[];
RUN(p:end,:)=[];
SA(p:end,:)=[];
SL(p:end,:)=[];
SR(p:end,:)=[];
TSTC(p:end,:)=[];
TSTI(p:end,:)=[];
TSTO(p:end,:)=[];
V(p:end,:)=[];

% get rid of the the initial 12 psi runs in the 2nd file which are
% duplicates of later data
p = 46007; % data point which we want to cut off from
p2 = 47243; % data point which we want to cut off to
AMBTMP(p:p2,:)=[];
ET(p:p2,:)=[];
FX(p:p2,:)=[];
FY(p:p2,:)=[];
FZ(p:p2,:)=[];
IA(p:p2,:)=[];
MX(p:p2,:)=[];
MZ(p:p2,:)=[];
N(p:p2,:)=[];
NFX(p:p2,:)=[];
NFY(p:p2,:)=[];
P(p:p2,:)=[];
RE(p:p2,:)=[];
RL(p:p2,:)=[];
RST(p:p2,:)=[];
RUN(p:p2,:)=[];
SA(p:p2,:)=[];
SL(p:p2,:)=[];
SR(p:p2,:)=[];
TSTC(p:p2,:)=[];
TSTI(p:p2,:)=[];
TSTO(p:p2,:)=[];
V(p:p2,:)=[];

%{
% plot V, P, SA, SL, IA vs ET to verify the data is good before procceding
figure;
plot(ET, V)
figure;
plot(ET, P)
figure;
plot(ET, SA)
figure;
plot(ET, IA)
figure;
plot(ET, FZ)
figure;
plot(ET, SL)
%}

%% Determine how to break up data in individual SL sweeps

% Each SL sweep starts at positive max SL, goes to negative max SL, and
% then back up to positive max SL, which makes the zero crossings hard to
% use for isolating SL sweeps.
% find the peaks of the SL channel to isolate SL sweeps
m = 1:length(SL);
[sl_peaks, z] = findpeaks(SL, MinPeakDistance=375);

%{
% plot raw slip angle data, plot the spline and zeros
figure;
subplot(5,1,1)
hold on
plot(m, SL)
plot(z,sl_peaks,'o')
hold off
xlabel('Point Count')
ylabel('Slip Ratio')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff
subplot(5,1,2)
hold on
plot(m, FZ)
plot(peak_locations,FZ(peak_locations),'o')
xlabel('Point Count')
ylabel('Normal Load')
subplot(5,1,3)
hold on
plot(m, IA)
xlabel('Point Count')
ylabel('IA')
subplot(5,1,4)
hold on
plot(m, P)
xlabel('Point Count')
ylabel('P')
subplot(5,1,5)
hold on
plot(m, V)
xlabel('Point Count')
ylabel('V')
%} 

figure;
hold on
plot(m, SL)
plot(z,sl_peaks,'o')
hold off
xlabel('Point Count')
ylabel('Slip Ratio')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff

%{
% Mz seems quite noisy, probably due to the tire not being round
figure;
subplot(2,1,1)
hold on
plot(m, MZ)
hold off
xlabel('Point Count')
ylabel('MZ')
subplot(2,1,2)
hold on
plot(m, FZ)
hold off
xlabel('Point Count')
ylabel('FZ')
%}

%% Capture individual SL sweep data:
% Don't want to remember any data from apprevious processing session:
clear fmdata;
% Subset the data into individual SL sweeps:
q = 0;
for n=1:(length(z)-1) % every 2 peaks represents a full SL sweep

    % capture data
    sl=SL(z(n):z(n+1)); % continuous n to n+1 represents going through an SL sweep
    sa=SA(z(n):z(n+1)); % constant
    fz=FZ(z(n):z(n+1)); % constant
    fy=FY(z(n):z(n+1)); % continuous
    fx=FX(z(n):z(n+1)); % continuous
    mz=MZ(z(n):z(n+1)); % continuous
    mx=MX(z(n):z(n+1)); % continuous
    rl=RL(z(n):z(n+1)); % continuous
    ia=IA(z(n):z(n+1)); % constant
    pressure=P(z(n):z(n+1)); % constant
    % Now we have collected the tire channels for each full slip ratio sweep.

    %% normalize the Fy & Fx channel to account for Fz noise
    % for each transient sweep, Fz is held constant, so we normalize Fy
    % by multiplying it by the ratio between the Fz mean and the isntantaneous Fz.
    fy = (fy * mean(fz)) ./ fz; % element-wise divide
    fx = (fx * mean(fz)) ./ fz; % element-wise divide
    
    %{
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
    %}

    %% Spline fitting the continuous data to subset it with increments of slip ratio
    % with some tighter tension:
    sp_fy=csaps(sl,fy,.999);
    sp_fx=csaps(sl,fx,.9999);
    sp_mz=csaps(sl,mz,.99);
    sp_mx=csaps(sl,mx,.99);
    sp_rl=csaps(sl,rl,.99);

    %% Check out Segment 9
    % Just out of curiosity, what kind of data are we dealing with?
    if isequal(n,23)
        fprintf("hit segment 9")

        figure;
        subplot(4,1,1)
        hold on
        plot(sl,fx,'.','color',[.5 .5 .5])
        fnplt(sp_fx,'b')
        title({['Fz= ' num2str(round(mean(fz))) ' N' ' SA= ' num2str(round(mean(sa))) ' deg'];['IA= ' num2str(round(mean(ia))) '°, ' num2str(round(mean(pressure))) ' psi']})
        xlabel('Slip Ratio')
        ylabel('Long. Force')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(fx) max(fx)],'color','k')
        legend('Test Data','Fitted Data')

        subplot(4,1,2)
        hold on
        plot(sl,fy,'.','color',[.5 .5 .5])
        fnplt(sp_fy,'b')
        xlabel('Slip Ratio')
        ylabel('Lateral Force')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(fy) max(fy)],'color','k')
        
        subplot(4,1,3)
        hold on
        plot(sl,mz,'.','color',[.5 .5 .5])
        fnplt(sp_mz,'b')
        xlabel('Slip Ratio')
        ylabel('Aligning Moment')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(mz) max(mz)],'color','k')

        subplot(4,1,4)
        hold on
        plot(sl,mx,'.','color',[.5 .5 .5])
        fnplt(sp_mx,'b')
        xlabel('Slip Ratio')
        ylabel('Overturning Moment') 
        line([min(sl) max(sl)],[0 0],'color','k') 
        line([0 0],[min(mz) max(mz)],'color','k')
    end
    
    %{
    %% fmdata is the organized, fitted data
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
%% Add Fz = 0N condition to data
% We know that is Fz = 0, there cannot be any Fy, Mz, Mx, because there is
% no load on the tire. By adding this known fact into our data, any models 
% we fit to this data will be accurate at extreme load transfer, or when
% the car is on two wheels.

incls = unique(round(fmdata(:,2)));
nincls = length(incls); % number of distinct inclination angles
slips = unique(round(fmdata(:,1)));
nslips = length(slips); % number of distinct slip angles
press = unique(round(fmdata(:,7)));
npress = length(press); % number of distinct pressures

% Create an array of our Fz = 0 condition:
% columns: SA, IA, Fz, Fy, Mz, Mx, P

% for each distinct P, for each distinct IA, concatenate array of slip
% angles
zero_load = [];
for n = 1:length(press)
    for q = 1:length(incls)
        zero_load = [zero_load; [slips, (zeros(length(slips), 1) + incls(q)), zeros(length(slips), 1), zeros(length(slips), 1), zeros(length(slips), 1), zeros(length(slips), 1), (zeros(length(slips), 1) + press(n))]];
    end
end

% concatenate array to existing data
fmdata = [fmdata; zero_load];

%% Sort Formatted Data Array 1st by Pressure, Camber, Slip, and finally Fz
% Use |sortrows| to preserve the array correspondence.
fmdata = sortrows(fmdata,[7,2,1,3])

%% Save organized tire data
save("../tire_data/d2704/d2704_7in_formatted_data.mat", "fmdata");

%% Data Sets (Slip, Load, Camber, Pressure):
% transpose these into column arrays
incls = incls';
slips = slips';
press = press';

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
fy0ip = reshape(fmdata0ip(:,4),nloads,nslips)';
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