% test and verify vehicle.m works correctly
clc; clear; clc;

zr25 = vehicle('zr25_data.xlsx');

zr25.create_simulink_parameters();

%testparam = Simulink.Parameter(5);
