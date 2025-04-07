clear;
clc;
% G=tf([34100, 39700], [1 0])
% Gd=c2d(G, .01, 'zoh')
% 
% % matrix

% K_p=[25320 8344.4 2483.7 1000 34132 21997
%     25320 8300 2483.7 6993.3 34132 34452
%     25320 8300 2483.7 8642 35054 46685
%     25320 8300 6836.2 71300 44181 50763
%     25320 8300 6836.2 71300 48881 50370];
% K_i=[11396 11663 5679.9 1000 29683 66594
%     11396 12000 5679.9 8782.2 29683 70850
%     11396 12000 5679.9 2193.8 34602 83942
%     11396 12000 46016 13884 42895 93304
%     11396 12000 46016 13884 49547 101260];
% v=[2 6 12 20 30 50];
% deg=[10 20 40 60 90];

K_p=[25320 8344.4 2483.7 1000 34132 21997 25320 8300 2483.7 6993.3 34132 34452 25320 8300 2483.7 8642 35054 46685 25320 8300 6836.2 71300 44181 50763 25320 8300 6836.2 71300 48881 50370];
K_i=[11396 11663 5679.9 1000 29683 66594 11396 12000 5679.9 8782.2 29683 70850 11396 12000 5679.9 2193.8 34602 83942 11396 12000 46016 13884 42895 93304 11396 12000 46016 13884 49547 101260];
v=[2 6 12 20 30 50 2 6 12 20 30 50 2 6 12 20 30 50 2 6 12 20 30 50 2 6 12 20 30 50 ];
deg=[10 10 10 10 10 10 20 20 20 20 20 20 40 40 40 40 40 40 60 60 60 60 60 60 90 90 90 90 90 90];

[Xq, Yq] = meshgrid([2:1:50], [10:2:90]);

Vq1=griddata(v, deg, K_p, Xq, Yq, 'cubic');
Vq2=griddata(v, deg, K_i, Xq, Yq, 'cubic');

% Plot the original data and the interpolated result
figure;
subplot(1,2,1);
scatter(v, deg, 100, K_i, 'filled'); % Original scattered data points
title('Original Irregularly Spaced Data');
colorbar;

subplot(1,2,2);
surf(Xq, Yq, Vq2); % Interpolated surface
title('Interpolated Data');
xlabel('v');
ylabel('deg');
zlabel('Ki');
colorbar;


figure;
subplot(1,2,1);
scatter(v, deg, 100, K_p, 'filled'); % Original scattered data points
title('Original Irregularly Spaced Data');
colorbar;

subplot(1,2,2);
surf(Xq, Yq, Vq1); % Interpolated surface
title('Interpolated Data');
xlabel('v');
ylabel('deg');
zlabel('Kp');
colorbar;
