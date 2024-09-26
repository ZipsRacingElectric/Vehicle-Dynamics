% Note: OptimumLap outputs data at varying sample rate (like wtf). It averages to
% 80 Hz sampling rate, which means the FFT is really only useful up to 40
% Hz. Actual measured data would make this much more valid. There are
% spikes in the data probably from track section discontinuities, and being
% a point mass sim, no over/undershoot of yaw-rate, or side slip angle
% effects are captured in this data. Only yaw-rates developed from
% cornering speeds and the curvieness of FSAE tracks is captured.

clear; clc;

A = readmatrix("ZRE23 OptimumLap.csv");
time = A(:,1);
yaw_rate = A(:,2);

% sample frequency
Fs = 1 / 0.0125;
% length of signal
L = length(time);

% convert from rpm to rad/s
yaw_rate = yaw_rate * (2*pi / 1) * (1 / 60);

% FFT on yaw rate
yaw_fft = fft(yaw_rate);
P2 = abs(yaw_fft/L);
% single ended FFT
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs/L*(0:(L/2));
plot(f, P1, "LineWidth", 3)
title("Frequency Spectrum of Yaw Rate, Optimum Lap, 80Hz Sample Rate")
xlabel("f (Hz)");
ylabel("|P1(f)|");