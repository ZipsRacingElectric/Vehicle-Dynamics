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
% ensure repository is part of MATLAB path
directory = fileparts(which(mfilename)); 
parentDirectory = fileparts(directory);
addpath(genpath(parentDirectory));

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
p = 92461; % data point which we want to cut off to
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

% plot V, P, SA, SL, IA vs ET to verify the data is good before procceding
figure;
subplot(6,1,1)
plot(ET,V,'.')
title('Velocity vs ET')
xlabel('ET (s)')
ylabel('Velocity (m/s)')

subplot(6,1,2)
plot(ET,P,'.')
title('Pressure vs ET')
xlabel('ET (s)')
ylabel('Pressue (PSI)')

subplot(6,1,3)
plot(ET,SA,'.')
title('Slip Angle vs ET')
xlabel('ET (s)')
ylabel('Slip Angle (m/s)')

subplot(6,1,4)
plot(ET,IA,'.')
title('Inclination Angle vs ET')
xlabel('ET (s)')
ylabel('Inclination Angle (deg)')

subplot(6,1,5)
plot(ET,FZ,'.')
title('Normal Force vs ET')
xlabel('ET (s)')
ylabel('Normal Force (N)')

subplot(6,1,6)
plot(ET,SL,'.')
title('Slip Ratio vs ET')
xlabel('ET (s)')
ylabel('Slip Ratio')


% Note: You will find that some tests, at the most IA and highest Fz, the
% slip angle sweep pretty much starts at zero, sweeps negative and back to 
% zero, with not much positive SL data. Need to check these sweeps and make
% sure that our csaps fit is acceptable to infer the positivie data

%% Determine how to break up data in individual SL sweeps

% Each SL sweep starts at positive max SL, goes to negative max SL, and
% then back up to positive max SL, which makes the zero crossings hard to
% use for isolating SL sweeps.
% instead, find the peaks of the SL channel to isolate SL sweeps
% these points will allow us to find the start / end indexes of the data
% for each sweep
m = 1:length(SL);
[sl_peaks, z] = findpeaks(SL, MinPeakDistance=375);

figure;
hold on
plot(m, SL)
plot(z,sl_peaks,'o')
hold off
xlabel('Point Count')
ylabel('Slip Ratio')
legend('Test Data','Computed Slip Points of Interest'),legend Boxoff

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


%% Capture individual SL sweep data:
% Don't want to remember any data from apprevious processing session:
clear fmdata;

q = 0;
for n=1:(length(z)-1) % every 2 peaks represents a full SL sweep

    % capture data
    % alse we are eliminating the peak data points. Because the SL sweepss
    % are discontinuous, we dont want to fit a curve to data that is part
    % of another sweep
    sl=SL((z(n)+1):(z(n+1)-1)); % continuous n to n+1 represents going through an SL sweep
    sa=SA((z(n)+1):(z(n+1)-1)); % constant
    fz=FZ((z(n)+1):(z(n+1)-1)); % constant
    fy=FY((z(n)+1):(z(n+1)-1)); % continuous
    fx=FX((z(n)+1):(z(n+1)-1)); % continuous
    mz=MZ((z(n)+1):(z(n+1)-1)); % continuous
    mx=MX((z(n)+1):(z(n+1)-1)); % continuous
    rl=RL((z(n)+1):(z(n+1)-1)); % continuous
    ia=IA((z(n)+1):(z(n+1)-1)); % constant
    pressure=P(z(n):(z(n+1)-1)); % constant
    % Now we have collected the tire channels for each full slip ratio sweep.

    %% Deal with SL sweeps missing data
    % some sweeps never got the positive side of SR. Need enough data to
    % ensure the spline will predict a resonable shape. We will simply copy
    % the raw data about the SR = 0 axis and mirror it. Not ideal because
    % the data isn't necessarily symmetric but good enough for now.
    % This is okay because this only happens at low load, high IA where Fx
    % is symmetrical and Mx, Mz seem constant
    if (max(sl) < 0.1 )
        % remove data with positive SL
        positive_indexes = find(sl > 0);
        sl(positive_indexes) = [];
        sa(positive_indexes) = [];
        fz(positive_indexes) = [];
        fy(positive_indexes) = [];
        fx(positive_indexes) = [];
        mz(positive_indexes) = [];
        mx(positive_indexes) = [];
        rl(positive_indexes) = [];
        ia(positive_indexes) = [];
        pressure(positive_indexes) = [];

        % mirror negative SL data and concatenate it
        sl = [sl; -sl];
        sa = [sa; sa];
        fz = [fz; fz];
        fy = [fy; fy];
        fx = [fx; -fx];
        mz = [mz; mz];
        mx = [mx; mx];
        rl = [rl; rl];
        ia = [ia; ia];
        pressure = [pressure; pressure];
    end


    %% normalize the Fy & Fx channel to account for Fz noise
    % Note: this was not helpful as Fz noise did not show in Fx/Fy channels
    %{
    % for each transient sweep, Fz is held constant, so we normalize Fy
    % by multiplying it by the ratio between the Fz mean and the isntantaneous Fz.
    fy = (fy * mean(fz)) ./ fz; % element-wise divide
    fx = (fx * mean(fz)) ./ fz; % element-wise divide
    %}

    %% Filtering Mz
    % TODO

    %% Filtering Mz
    % TODO

    %% Spline fitting the continuous data to subset it with increments of slip ratio
    % the smoothing parameters are adjusted to get a quality fit
    % without overfitting outliers
    sp_fy = csaps(sl,fy,.9999);
    sp_fx = csaps(sl,fx,.9999);
    sp_mz = csaps(sl,mz,.9999);
    sp_mx = csaps(sl,mx,.9999);
    sp_rl = csaps(sl,rl,.99);

    %% Check out Segment 9
    % Just out of curiosity, what kind of data are we dealing with?
    if isequal(n,95)

        figure;
        subplot(4,1,1)
        hold on
        plot(sl,fx,'.','color',[.5 .5 .5])
        fnplt(sp_fx,'b')
        title({['Fz= ' num2str(round(mean(fz))) ' N' ' SA= ' num2str(round(mean(sa))) ' deg'];['IA= ' num2str(round(mean(ia))) '°, ' num2str(round(mean(pressure))) ' psi']})
        xlabel('Slip Ratio')
        ylabel('Fx (N)')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(fx) max(fx)],'color','k')
        legend('Test Data','Fitted Data')
        % set(gca,'XTick',[],'YTick',[]) % for public repository / reports

        subplot(4,1,2)
        hold on
        plot(sl,fy,'.','color',[.5 .5 .5])
        fnplt(sp_fy,'b')
        xlabel('Slip Ratio')
        ylabel('Fy (N)')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(fy) max(fy)],'color','k')
        % set(gca,'XTick',[],'YTick',[])
        
        subplot(4,1,3)
        hold on
        plot(sl,mz,'.','color',[.5 .5 .5])
        fnplt(sp_mz,'b')
        xlabel('Slip Ratio')
        ylabel('Aligning Moment Mz (N-m)')
        line([min(sl) max(sl)],[0 0],'color','k')
        line([0 0],[min(mz) max(mz)],'color','k')
        % set(gca,'XTick',[],'YTick',[])

        subplot(4,1,4)
        hold on
        plot(sl,mx,'.','color',[.5 .5 .5])
        fnplt(sp_mx,'b')
        xlabel('Slip Ratio')
        ylabel('Overturning Moment Mx (N-m)') 
        line([min(sl) max(sl)],[0 0],'color','k') 
        line([0 0],[min(mz) max(mz)],'color','k')
        % set(gca,'XTick',[],'YTick',[])
    end
    
    %% fmdata is the organized, fitted data
    % columns: SL, SA IA, P, Fz, Fy, Mz, Mx

    step_size = 0.02;
    start_slip = floor(min(sl) / step_size) * step_size; % Smallest number <= min(sl) evenly divisible by step size
    end_slip = ceil(max(sl) / step_size) * step_size; % Smallest number >= max(sl) evenly divisible by step size
    for slip = start_slip:step_size:end_slip % evaluate data and fitted splines for slip ratios in small increments
        q = q+1;
        fmdata(q,1) = slip;
        fmdata(q,2) = round(mean(sa)); 
        fmdata(q,3) = round(mean(ia));
        fmdata(q,4) = round(mean(pressure));
        fmdata(q,5) = mean(fz);
        fmdata(q,6) = fnval(sp_fy,slip);
        fmdata(q,7) = fnval(sp_mz,slip);
        fmdata(q,8) = fnval(sp_mx,slip);
    end 
end

%% Add Fz = 0N condition to data
% We know that at Fz = 0, there cannot be any Fx, Fy, Mz, Mx, because there is
% no load on the tire. By adding this known fact into our data, any models 
% we fit to this data will be accurate at extreme load transfer, or when
% the car is on two wheels.

srs = unique(round(fmdata(:,1), 2));
nsrs = length(srs); % number of distinct slip ratios
sas = unique(round(fmdata(:,2)));
nsas = length(sas); % number of distinct slip angles
incls = unique(round(fmdata(:,3)));
nincls = length(incls); % number of distinct inclination angles
press = unique(round(fmdata(:,4)));
npress = length(press); % number of distinct pressures

% Create an array of our Fz = 0 condition:
% columns: SL, SA, IA, P, Fz, Fy, Mz, Mx

% for each distinct P, IA, and SA concatenate array of slip ratios
zero_load = [];
for n = 1:length(press)
    for q = 1:length(incls)
        for i = 1:length(sas)
            zero_load = [zero_load; [srs, (zeros(length(srs), 1) + sas(q)), (zeros(length(srs), 1) + incls(q)), (zeros(length(srs), 1) + press(n)), zeros(length(srs), 1), zeros(length(srs), 1), zeros(length(srs), 1), zeros(length(srs), 1)]];
        end
    end
end

% concatenate array to existing data
fmdata = [fmdata; zero_load];

% remove any redundant rows in the data
fmdata = unique(fmdata, 'rows');

%% Sort Formatted Data Array by Pressure, Camber, Slip Angle, Slip Ratio, and finally Fz
% Use |sortrows| to preserve the array correspondence.
fmdata = sortrows(fmdata,[4,3,2,1, 5])

%% Save organized tire data
save("../tire_data/d2704/d2704_7in_formatted_combined_data.mat", "fmdata");
