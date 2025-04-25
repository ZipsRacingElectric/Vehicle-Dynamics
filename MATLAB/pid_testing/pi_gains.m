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
surf(Xq, Yq, Vq2); % Interpolated surface
title('Ki Gain Map');
xlabel('Velocity (m/s)');
ylabel('Steering Angle (degrees)');
zlabel('Ki');
colorbar;

figure;
surf(Xq, Yq, Vq1); % Interpolated surface
title('Kp Gain Map');
xlabel('Velocity (m/s)');
ylabel('Steering Angle (degrees)');
zlabel('Kp');
colorbar;

% Improved Ki Gain Map
figure;
surf(Xq, Yq, Vq2);                    % Smooth the color transitions
colormap(parula);                     % Use a visually pleasing colormap
colorbar;
title('Ki Gain Map', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Velocity (m/s)', 'FontSize', 12);
ylabel('Steering Angle (degrees)', 'FontSize', 12);
zlabel('Ki', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);   
colorbar;% Set a better viewing angle

% Improved Kp Gain Map
figure;
surf(Xq, Yq, Vq1);
colormap(parula);                      % A slightly different colormap
colorbar;
title('Kp Gain Map', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Velocity (m/s)', 'FontSize', 12);
ylabel('Steering Angle (degrees)', 'FontSize', 12);
zlabel('Kp', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% Save pi gains to .mat file
save('results/tuned_gains.mat', K_p, K_i, v, deg)
