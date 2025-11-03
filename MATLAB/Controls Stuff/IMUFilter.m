% IMU data filtering test

clc, clear

data = readtable('C:\Users\ATuck\Documents\GitHub\Vehicle-Dynamics\MATLAB\abbieTest\IMUDataForFiltering.csv');

fs = 100;               % Sampling rate in Hz
fc = 10;                  % Desired cutoff frequency in Hz (adjust as needed)
order = 4;              % Filter order (4th-order Butterworth)

% define data variables 
IMU_X = data.CAN1_ECUMASTER_GPS_IMU1_X_ACCELERATION ;
IMU_Y = data.CAN1_ECUMASTER_GPS_IMU1_Y_ACCELERATION ;
IMU_Z = data.CAN1_ECUMASTER_GPS_IMU1_Z_ACCELERATION ;
SAS = data.CAN1_VCU_SENSOR_INPUT_STEERING_ANGLE ;
longitude = data.CAN1_ECUMASTER_GPS_POSITION_LONGITUDE ;
Latitude = data.CAN1_ECUMASTER_GPS_POSITION_LATITUDE ;
time = data.timestamps ;

%% ------------------- Design Low-Pass Filter -------------------
[b,a] = butter(order, fc/(fs/2), 'low');   % Butterworth low-pass

%% ------------------- Apply Zero-Phase Filtering -------------------
IMUX_f = filtfilt(b, a, IMU_X);
IMUY_f = filtfilt(b, a, IMU_Y);
IMUZ_f = filtfilt(b, a, IMU_Z);


figure(1)
subplot(3,1,1) 
plot(time(9000:10000, 1), IMUX_f(9000:10000, 1), Marker='🐸')
subplot(3,1,2)
%plot3(longitude(8000:10000,1),Latitude(8000:10000,1),time(8000:10000,1))
grid on


%figure(2)
%plot3(IMUX_f,IMUY_f, IMUZ_f)
%grid on